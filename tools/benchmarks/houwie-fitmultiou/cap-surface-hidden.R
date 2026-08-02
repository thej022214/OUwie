#!/usr/bin/env Rscript

# The near-optimum accuracy of the capped likelihood for the two model classes
# the cap-by-model study left open: four-regime hidden-state models, and the
# fully general OUMVA/OUMA parameterizations.
#
# The reference here is the pruning likelihood at max_components = Inf rather
# than exhaustive enumeration. That is exact for the stated resolution, and the
# convergence study already checked it against enumeration to ~1e-14 for these
# same model classes; using it as the reference is what makes a scan of this many
# parameter points affordable at four regimes.
#
# The statistic that matters is the error *near the optimum*. A likelihood that
# is hundreds of log units wrong six hundred log units below the optimum costs
# nothing, because no optimizer goes there and the error is negative, which makes
# an already bad region look worse. What costs something is an error that moves
# the optimum, or a positive error that manufactures a mode.

benchmark_dir <- normalizePath(dirname(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
))
repo_root <- normalizePath(file.path(benchmark_dir, "..", "..", ".."))
suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE,
                    helpers = FALSE)
})
source(file.path(benchmark_dir, "mixture-pruning-prototype.R"))

out_dir <- file.path(benchmark_dir, "results", "cap-surface")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
caps <- c(4, 8, 16, 32)

make_case <- function(n_tip, seed, q = 0.3){
  set.seed(seed)
  phy <- pbtree(n = n_tip, scale = 10)
  Q <- matrix(c(-q, q, q, -q), 2, 2, byrow = TRUE)
  sim <- hOUwie.sim(phy, Q = Q, root.freqs = c(1, 0), alpha = c(0.5, 0.5),
                    sigma.sq = c(0.2, 0.2), theta0 = 5, theta = c(5, 10))
  phy <- reorder.phylo(phy, "pruningwise")
  values <- setNames(sim$data[, 3], sim$data[, 1])[phy$tip.label]
  observed <- setNames(as.integer(sim$data[, 2]), sim$data[, 1])[phy$tip.label]
  allowed <- lapply(observed, function(s) s)
  names(allowed) <- names(observed)
  list(phy = phy, values = values, allowed = allowed, observed = observed)
}

# Every tip is ambiguous between the two rate categories of its observed state,
# which is what makes a hidden-state message carry four live components per node
# instead of one.
hide <- function(case){
  allowed <- lapply(case$observed,
                    function(s) if(s == 1L) c(1L, 2L) else c(3L, 4L))
  names(allowed) <- names(case$observed)
  list(phy = case$phy, values = case$values, allowed = allowed)
}

Q4 <- matrix(0, 4, 4)
Q4[1, 2] <- Q4[2, 1] <- Q4[3, 4] <- Q4[4, 3] <- 0.15
Q4[1, 3] <- Q4[3, 1] <- Q4[2, 4] <- Q4[4, 2] <- 0.30
diag(Q4) <- -rowSums(Q4)
Q2 <- matrix(c(-0.3, 0.3, 0.3, -0.3), 2, 2, byrow = TRUE)

studies <- list(
  # four regimes, tips ambiguous between rate categories
  list(label = "hidden theta separation", tips = 8, Q = Q4, root = rep(.25, 4),
       hidden = TRUE, axis = seq(0, 14, length.out = 25),
       build = function(v) list(alpha = c(.25, 1, .25, 1),
                                sigma = c(.15, .6, .15, .6),
                                theta = c(5, 5, 5 + v, 5 + v))),
  list(label = "hidden alpha ratio", tips = 8, Q = Q4, root = rep(.25, 4),
       hidden = TRUE, axis = exp(seq(log(1), log(20), length.out = 25)),
       build = function(v) list(alpha = c(.25, .25 * v, .25, .25 * v),
                                sigma = c(.3, .3, .3, .3),
                                theta = c(5, 5, 12, 12))),
  list(label = "hidden sigma.sq ratio", tips = 8, Q = Q4, root = rep(.25, 4),
       hidden = TRUE, axis = exp(seq(log(1), log(20), length.out = 25)),
       build = function(v) list(alpha = c(.5, .5, .5, .5),
                                sigma = c(.15, .15 * v, .15, .15 * v),
                                theta = c(5, 5, 12, 12))),
  # OUMVA and OUMA: everything free at once, which is the case the cap-by-model
  # study flagged and the recovery suite never touched
  list(label = "OUMVA joint", tips = 12, Q = Q2, root = c(.5, .5),
       hidden = FALSE, axis = exp(seq(log(1), log(16), length.out = 25)),
       build = function(v) list(alpha = c(.25, .25 * v),
                                sigma = c(.15, .15 * v),
                                theta = c(5, 5 + 12 * (v - 1) / 15))),
  list(label = "OUMA joint", tips = 12, Q = Q2, root = c(.5, .5),
       hidden = FALSE, axis = exp(seq(log(1), log(16), length.out = 25)),
       build = function(v) list(alpha = c(.25, .25 * v), sigma = c(.3, .3),
                                theta = c(5, 5 + 12 * (v - 1) / 15))),
  list(label = "OUMVA theta separation", tips = 12, Q = Q2, root = c(.5, .5),
       hidden = FALSE, axis = seq(0, 14, length.out = 25),
       build = function(v) list(alpha = c(.25, 1), sigma = c(.15, .6),
                                theta = c(5, 5 + v)))
)

rows <- list()
for(study in studies){
  case <- make_case(study$tips, 5200 + study$tips)
  if(study$hidden) case <- hide(case)
  cat(sprintf("%s (tips=%d, K=%d)\n", study$label, study$tips, nrow(study$Q)))
  evaluate <- function(cap){
    vapply(study$axis, function(v){
      p <- study$build(v)
      mixture_pruning_loglik(case$phy, case$allowed, case$values, study$Q,
                             p$alpha, p$sigma, p$theta, root_p = study$root,
                             resolution = 1L, max_components = cap)
    }, numeric(1))
  }
  began <- Sys.time()
  exact <- evaluate(Inf)
  cat(sprintf("  reference (Inf) took %.1fs\n",
              as.numeric(difftime(Sys.time(), began, units = "secs"))))
  utils::flush.console()
  for(cap in caps){
    approx <- evaluate(cap)
    err <- approx - exact
    below <- max(exact) - exact
    near <- function(limit) if(any(below <= limit)){
      max(abs(err[below <= limit]))
    }else NA_real_
    modes <- function(y){
      i <- seq_along(y)[-c(1, length(y))]
      i[y[i] > y[i - 1] & y[i] > y[i + 1]]
    }
    same_modes <- identical(modes(exact), modes(approx))
    rows[[length(rows) + 1L]] <- data.frame(
      label = study$label, tips = study$tips, regimes = nrow(study$Q),
      cap = cap, err_within2 = near(2), err_within10 = near(10),
      err_all = max(abs(err)), max_positive = max(err),
      argmax_exact = study$axis[which.max(exact)],
      argmax_approx = study$axis[which.max(approx)],
      same_mode_structure = same_modes, stringsAsFactors = FALSE)
    cat(sprintf("  cap=%-3d within2=%8.2e within10=%8.2e all=%9.3f  modes %s\n",
                cap, near(2), near(10), max(abs(err)),
                if(same_modes) "match" else "DIFFER"))
    utils::flush.console()
  }
}
rows <- do.call(rbind, rows)
write.csv(rows, file.path(out_dir, "hidden-oumva-scans.csv"), row.names = FALSE)
cat("\nwrote", file.path(out_dir, "hidden-oumva-scans.csv"), "\n")
