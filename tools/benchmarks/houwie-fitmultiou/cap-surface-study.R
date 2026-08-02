#!/usr/bin/env Rscript

# Is the capped likelihood usable as an objective function?
#
# The convergence study varied the cap at a fixed parameter vector and asked how
# large the error was. That is the wrong question for a method whose whole job is
# to be optimized. An objective can carry a large error and still find the right
# optimum, provided the error is a smooth function of the parameters; and it can
# carry a small error and be useless, if the error jumps as components are merged
# differently either side of a step.
#
# So here the cap is fixed and the parameters are scanned, against exhaustive
# enumeration of the same midpoint model. Reported per scan:
#
#   max |error|          how wrong the surface is
#   max |d error|        the largest jump between adjacent grid points, which is
#                        what a gradient-free optimizer walks into
#   argmax displacement  where the capped surface puts its optimum along the scan
#                        relative to where the exact surface puts it
#
# Trees are kept to twelve tips so that max_components = Inf is exact and cheap.

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

plan_for <- function(phy, n_states){
  OUwie:::getHOUwieTreePlan(
    phy,
    OUwie:::getEdgeLiks(phy,
                        data.frame(sp = phy$tip.label, reg = 1),
                        n_states, 1, max(branching.times(phy)) + 1)
  )
}

enumerate_midpoint <- function(phy, tip_allowed, tip_values, Q, alpha,
                               sigma_sq, theta, root_p){
  n_tip <- Ntip(phy)
  data_ou <- data.frame(sp = phy$tip.label, reg = 1,
                        x = unname(tip_values[phy$tip.label]),
                        stringsAsFactors = FALSE)
  tree_plan <- plan_for(phy, nrow(Q))
  transitions <- lapply(phy$edge.length, function(t) expm::expm(Q * t))
  grid <- expand.grid(c(tip_allowed[phy$tip.label],
                        rep(list(seq_len(nrow(Q))), Nnode(phy))))
  terms <- numeric(nrow(grid))
  for(i in seq_len(nrow(grid))){
    states <- as.integer(grid[i, ])
    discrete <- log(root_p[states[n_tip + 1L]]) +
      sum(vapply(seq_len(nrow(phy$edge)), function(e){
        log(transitions[[e]][states[phy$edge[e, 1]], states[phy$edge[e, 2]]])
      }, numeric(1)))
    midpoint <- mapply(function(a, d, len){
      setNames(rep(len / 2, 2), c(states[a], states[d]))
    }, a = phy$edge[, 1], d = phy$edge[, 2], len = phy$edge.length,
    SIMPLIFY = FALSE)
    terms[i] <- discrete + OUwie:::OUwie.basic(
      phy, data_ou, simmap.tree = TRUE, scaleHeight = FALSE, alpha = alpha,
      sigma.sq = sigma_sq, theta = theta, algorithm = "three.point",
      tip.paths = tree_plan$tip.paths, tip.fog = "none", map = midpoint,
      tree.plan = tree_plan, map.states = as.character(seq_len(nrow(Q))))
  }
  log_sum_exp(terms)
}

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
  list(phy = phy, values = values, allowed = allowed)
}

Q2 <- matrix(c(-0.3, 0.3, 0.3, -0.3), 2, 2, byrow = TRUE)
caps <- c(4, 8, 16, 32)

# Each scan varies one parameter along a grid, holding the rest fixed. The
# builder returns the (alpha, sigma.sq, theta) triple for a scan point.
scans <- list(
  list(name = "theta separation", axis = seq(0, 14, length.out = 29),
       build = function(v) list(alpha = c(0.5, 0.5), sigma = c(0.3, 0.3),
                                theta = c(5, 5 + v))),
  list(name = "alpha (shared)", axis = exp(seq(log(0.02), log(4),
                                               length.out = 29)),
       build = function(v) list(alpha = c(v, v), sigma = c(0.3, 0.3),
                                theta = c(5, 12))),
  list(name = "alpha ratio", axis = exp(seq(log(1), log(20),
                                            length.out = 29)),
       build = function(v) list(alpha = c(0.25, 0.25 * v),
                                sigma = c(0.3, 0.3), theta = c(5, 12))),
  list(name = "sigma.sq ratio", axis = exp(seq(log(1), log(20),
                                               length.out = 29)),
       build = function(v) list(alpha = c(0.5, 0.5),
                                sigma = c(0.15, 0.15 * v),
                                theta = c(5, 12)))
)

rows <- list()
for(n_tip in c(8, 12)){
  case <- make_case(n_tip, 5200 + n_tip)
  for(scan in scans){
    cat(sprintf("tips=%d  scan=%s\n", n_tip, scan$name))
    exact <- vapply(scan$axis, function(v){
      p <- scan$build(v)
      enumerate_midpoint(case$phy, case$allowed, case$values, Q2, p$alpha,
                         p$sigma, p$theta, c(0.5, 0.5))
    }, numeric(1))
    for(cap in caps){
      approx <- vapply(scan$axis, function(v){
        p <- scan$build(v)
        mixture_pruning_loglik(case$phy, case$allowed, case$values, Q2,
                               p$alpha, p$sigma, p$theta, root_p = c(0.5, 0.5),
                               resolution = 1L, max_components = cap)
      }, numeric(1))
      rows[[length(rows) + 1L]] <- data.frame(
        tips = n_tip, scan = scan$name, cap = cap,
        point = seq_along(scan$axis), value = scan$axis,
        exact = exact, approx = approx, error = approx - exact,
        stringsAsFactors = FALSE)
      err <- approx - exact
      cat(sprintf("  cap=%-3d max|err|=%9.3f  max|jump in err|=%9.3f  "  ,
                  cap, max(abs(err)), max(abs(diff(err)))))
      cat(sprintf("argmax exact=%.3g approx=%.3g\n",
                  scan$axis[which.max(exact)], scan$axis[which.max(approx)]))
      utils::flush.console()
    }
  }
}
rows <- do.call(rbind, rows)
write.csv(rows, file.path(out_dir, "parameter-scans.csv"), row.names = FALSE)

summary_rows <- do.call(rbind, lapply(split(rows, list(rows$tips, rows$scan,
                                                       rows$cap), drop = TRUE),
                                      function(d){
  d <- d[order(d$point), ]
  data.frame(tips = d$tips[1], scan = d$scan[1], cap = d$cap[1],
             max_abs_error = max(abs(d$error)),
             max_error_jump = max(abs(diff(d$error))),
             argmax_exact = d$value[which.max(d$exact)],
             argmax_approx = d$value[which.max(d$approx)],
             stringsAsFactors = FALSE)
}))
write.csv(summary_rows, file.path(out_dir, "scan-summary.csv"),
          row.names = FALSE)
cat("\nwrote", file.path(out_dir, "scan-summary.csv"), "\n")
