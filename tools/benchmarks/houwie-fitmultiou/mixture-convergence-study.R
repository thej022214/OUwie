#!/usr/bin/env Rscript

# Choose the mixture-reduction defaults from evidence rather than from a guess.
#
# Every study below measures error against exhaustive enumeration of the same
# midpoint history model, so the reference is exact. The decision-relevant
# question is not "is cap 4 accurate on six tips" but "does the error compound
# with depth", because the trees this has to run on have three hundred tips.
#
#   1. cap x model class       - which parameterizations are hard to merge
#   2. cap x tree size         - how the error grows with depth
#   3. cap x regime count      - whether the cap has to scale with K
#   4. resolution correctness  - validate the resolution > 1 code path, which
#                                nothing has checked yet
#   5. resolution convergence  - self-convergence on the 300-tip fixture
#
# Results are written after each study so a long run can be inspected early.

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

out_dir <- file.path(benchmark_dir, "results", "mixture-convergence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
caps <- c(2, 4, 8, 16, 32, Inf)

plan_for <- function(phy, n_states){
  OUwie:::getHOUwieTreePlan(
    phy,
    OUwie:::getEdgeLiks(phy,
                        data.frame(sp = phy$tip.label, reg = 1),
                        n_states, 1, max(branching.times(phy)) + 1)
  )
}

# ---------------------------------------------------------------------------
# Exhaustive enumeration of the midpoint model: one regime point per node, the
# edge painted half with the ancestor's regime and half with the descendant's.
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# The same enumeration at resolution 2: an extra regime point in the middle of
# every edge, so each edge carries quarter/half/quarter painting and two
# transition steps. This exists to validate the resolution > 1 propagation.
# ---------------------------------------------------------------------------
enumerate_resolution_two <- function(phy, tip_allowed, tip_values, Q, alpha,
                                     sigma_sq, theta, root_p){
  n_tip <- Ntip(phy)
  n_edge <- nrow(phy$edge)
  data_ou <- data.frame(sp = phy$tip.label, reg = 1,
                        x = unname(tip_values[phy$tip.label]),
                        stringsAsFactors = FALSE)
  tree_plan <- plan_for(phy, nrow(Q))
  half_steps <- lapply(phy$edge.length,
                       function(t) expm::expm(Q * (t / 2)))
  node_grid <- c(tip_allowed[phy$tip.label],
                 rep(list(seq_len(nrow(Q))), Nnode(phy)))
  grid <- expand.grid(c(node_grid, rep(list(seq_len(nrow(Q))), n_edge)))
  terms <- numeric(nrow(grid))
  for(i in seq_len(nrow(grid))){
    row <- as.integer(grid[i, ])
    states <- row[seq_len(n_tip + Nnode(phy))]
    interior <- row[(n_tip + Nnode(phy)) + seq_len(n_edge)]
    discrete <- log(root_p[states[n_tip + 1L]])
    painted <- vector("list", n_edge)
    for(e in seq_len(n_edge)){
      from <- states[phy$edge[e, 1]]
      mid <- interior[e]
      to <- states[phy$edge[e, 2]]
      discrete <- discrete + log(half_steps[[e]][from, mid]) +
        log(half_steps[[e]][mid, to])
      quarter <- phy$edge.length[e] / 4
      painted[[e]] <- setNames(c(quarter, 2 * quarter, quarter),
                               c(from, mid, to))
    }
    terms[i] <- discrete + OUwie:::OUwie.basic(
      phy, data_ou, simmap.tree = TRUE, scaleHeight = FALSE, alpha = alpha,
      sigma.sq = sigma_sq, theta = theta, algorithm = "three.point",
      tip.paths = tree_plan$tip.paths, tip.fog = "none", map = painted,
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

sweep_caps <- function(label, case, Q, alpha, sigma_sq, theta, root_p,
                       enumerator = enumerate_midpoint, resolution = 1L,
                       extra = list()){
  exact <- enumerator(case$phy, case$allowed, case$values, Q, alpha, sigma_sq,
                      theta, root_p)
  do.call(rbind, lapply(caps, function(cap){
    began <- Sys.time()
    pruned <- mixture_pruning_loglik(case$phy, case$allowed, case$values, Q,
                                     alpha, sigma_sq, theta, root_p = root_p,
                                     resolution = resolution,
                                     max_components = cap)
    row <- data.frame(setting = label, tips = Ntip(case$phy),
                      regimes = nrow(Q), resolution = resolution,
                      cap = cap, enumerated = exact, pruned = pruned,
                      error = pruned - exact,
                      seconds = as.numeric(difftime(Sys.time(), began,
                                                    units = "secs")),
                      stringsAsFactors = FALSE)
    cbind(row, as.data.frame(extra, stringsAsFactors = FALSE))
  }))
}

announce <- function(rows){
  for(i in seq_len(nrow(rows))){
    cat(sprintf("  %-26s tips=%2d K=%d cap=%-4s error=%11.3e  %6.2fs\n",
                rows$setting[i], rows$tips[i], rows$regimes[i],
                ifelse(is.finite(rows$cap[i]), rows$cap[i], "Inf"),
                rows$error[i], rows$seconds[i]))
  }
  utils::flush.console()
}

# ---------------------------------------------------------------------------
# Study 1: cap x model class. Contrast scales whichever parameters the model
# lets vary, so "high" means well separated optima for OUM and a large rate
# ratio for OUV and OUA.
# ---------------------------------------------------------------------------
cat("study 1: cap x model class (8 tips, 2 regimes)\n")
case8 <- make_case(8, 5101)
contrasts <- list(low = 1, medium = 2, high = 3)
model_grid <- list(
  OUM   = function(s) list(alpha = c(0.5, 0.5), sigma = c(0.3, 0.3),
                           theta = c(5, 5 + c(1, 4, 12)[s])),
  OUV   = function(s) list(alpha = c(0.5, 0.5),
                           sigma = c(0.15, 0.15 * c(1.5, 4, 16)[s]),
                           theta = c(7, 7)),
  OUA   = function(s) list(alpha = c(0.25, 0.25 * c(1.5, 4, 16)[s]),
                           sigma = c(0.3, 0.3), theta = c(7, 7)),
  OUMV  = function(s) list(alpha = c(0.5, 0.5),
                           sigma = c(0.15, 0.15 * c(1.5, 4, 16)[s]),
                           theta = c(5, 5 + c(1, 4, 12)[s])),
  OUMA  = function(s) list(alpha = c(0.25, 0.25 * c(1.5, 4, 16)[s]),
                           sigma = c(0.3, 0.3),
                           theta = c(5, 5 + c(1, 4, 12)[s])),
  OUMVA = function(s) list(alpha = c(0.25, 0.25 * c(1.5, 4, 16)[s]),
                           sigma = c(0.15, 0.15 * c(1.5, 4, 16)[s]),
                           theta = c(5, 5 + c(1, 4, 12)[s])),
  BMS   = function(s) list(alpha = c(1e-10, 1e-10),
                           sigma = c(0.15, 0.15 * c(1.5, 4, 16)[s]),
                           theta = c(5, 5))
)
Q2 <- matrix(c(-0.3, 0.3, 0.3, -0.3), 2, 2, byrow = TRUE)
study1 <- list()
for(model in names(model_grid)){
  for(level in names(contrasts)){
    pars <- model_grid[[model]](contrasts[[level]])
    rows <- sweep_caps(paste(model, level), case8, Q2, pars$alpha, pars$sigma,
                       pars$theta, c(0.5, 0.5),
                       extra = list(model = model, contrast = level))
    announce(rows)
    study1[[length(study1) + 1L]] <- rows
  }
}
study1 <- do.call(rbind, study1)
write.csv(study1, file.path(out_dir, "cap-by-model.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Study 2: cap x tree size. This is the axis that decides whether a default is
# safe at three hundred tips.
# ---------------------------------------------------------------------------
cat("study 2: cap x tree size (2 regimes, high contrast)\n")
study2 <- list()
for(n_tip in c(6, 8, 10, 12)){
  case <- make_case(n_tip, 5200 + n_tip)
  for(model in c("OUM", "OUA")){
    pars <- model_grid[[model]](3)
    rows <- sweep_caps(paste(model, "high", n_tip), case, Q2, pars$alpha,
                       pars$sigma, pars$theta, c(0.5, 0.5),
                       extra = list(model = model, contrast = "high"))
    announce(rows)
    study2[[length(study2) + 1L]] <- rows
  }
}
study2 <- do.call(rbind, study2)
write.csv(study2, file.path(out_dir, "cap-by-tree-size.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Study 3: cap x regime count, using a hidden-state model where every tip is
# ambiguous between the two rate categories of its observed state.
# ---------------------------------------------------------------------------
cat("study 3: cap x regime count (5 tips, 4 regimes)\n")
case5 <- make_case(5, 5301)
observed5 <- vapply(case5$allowed, function(x) x[1], integer(1))
allowed4 <- lapply(observed5, function(s) if(s == 1L) c(1L, 2L) else c(3L, 4L))
names(allowed4) <- names(case5$allowed)
case5_hidden <- list(phy = case5$phy, values = case5$values,
                     allowed = allowed4)
Q4 <- matrix(0, 4, 4)
Q4[1, 2] <- Q4[2, 1] <- Q4[3, 4] <- Q4[4, 3] <- 0.15
Q4[1, 3] <- Q4[3, 1] <- Q4[2, 4] <- Q4[4, 2] <- 0.30
diag(Q4) <- -rowSums(Q4)
study3 <- list()
for(level in names(contrasts)){
  s <- contrasts[[level]]
  separation <- c(1, 4, 12)[s]
  ratio <- c(1.5, 4, 16)[s]
  rows <- sweep_caps(paste("hidden", level), case5_hidden, Q4,
                     alpha = c(0.25, 0.25 * ratio, 0.25, 0.25 * ratio),
                     sigma_sq = c(0.15, 0.15 * ratio, 0.15, 0.15 * ratio),
                     theta = c(5, 5, 5 + separation, 5 + separation),
                     root_p = rep(0.25, 4),
                     extra = list(model = "hidden", contrast = level))
  announce(rows)
  study3[[length(study3) + 1L]] <- rows
}
study3 <- do.call(rbind, study3)
write.csv(study3, file.path(out_dir, "cap-by-regime-count.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Study 4: is the resolution > 1 propagation correct? Enumerate the resolution
# 2 model directly and compare. Nothing has checked this path before.
# ---------------------------------------------------------------------------
cat("study 4: resolution 2 correctness (5 tips, 2 regimes)\n")
case5b <- make_case(5, 5401)
study4 <- list()
for(model in c("OUM", "OUMVA")){
  pars <- model_grid[[model]](3)
  rows <- sweep_caps(paste("res2", model), case5b, Q2, pars$alpha, pars$sigma,
                     pars$theta, c(0.5, 0.5),
                     enumerator = enumerate_resolution_two, resolution = 2L,
                     extra = list(model = model, contrast = "high"))
  announce(rows)
  study4[[length(study4) + 1L]] <- rows
}
study4 <- do.call(rbind, study4)
write.csv(study4, file.path(out_dir, "resolution-correctness.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# Study 5: resolution self-convergence on the fixture that started all this.
# There is no enumeration oracle here, so the question is whether the value
# stops moving.
# ---------------------------------------------------------------------------
cat("study 5: resolution self-convergence (300 tips)\n")
fixture <- readRDS(file.path(benchmark_dir, "results", "scaling",
                             "fixture-tips300-seed8.rds"))
phy300 <- reorder.phylo(fixture$phy, "pruningwise")
states300 <- as.integer(fixture$discrete)
names(states300) <- names(fixture$discrete)
study5 <- list()
for(resolution in c(1L, 2L, 4L, 8L, 16L)){
  began <- Sys.time()
  value <- mixture_pruning_loglik(phy300, states300, fixture$trait, Q2,
                                  alpha = rep(0.5, 2), sigma_sq = rep(0.2, 2),
                                  theta = c(5, 10), root_p = c(0.5, 0.5),
                                  resolution = resolution, max_components = 8L)
  elapsed <- as.numeric(difftime(Sys.time(), began, units = "secs"))
  cat(sprintf("  resolution=%2d  lnL=%12.4f  %6.1fs\n", resolution, value,
              elapsed))
  utils::flush.console()
  study5[[length(study5) + 1L]] <- data.frame(resolution = resolution,
                                              cap = 8, loglik = value,
                                              seconds = elapsed)
}
study5 <- do.call(rbind, study5)
write.csv(study5, file.path(out_dir, "resolution-convergence-300.csv"),
          row.names = FALSE)

cat("\nsmallest cap reaching each error tolerance, by setting:\n")
summary_rows <- do.call(rbind, lapply(
  split(rbind(study1, study2, study3), ~ setting), function(block){
    block <- block[order(block$cap), ]
    smallest <- function(tolerance){
      hit <- which(abs(block$error) < tolerance & is.finite(block$cap))
      if(!length(hit)) NA_real_ else block$cap[hit[1]]
    }
    data.frame(setting = block$setting[1], tips = block$tips[1],
               regimes = block$regimes[1],
               cap_for_0.1 = smallest(0.1), cap_for_0.01 = smallest(0.01),
               error_at_8 = block$error[block$cap == 8],
               stringsAsFactors = FALSE)
  }))
print(summary_rows, row.names = FALSE, digits = 3)
write.csv(summary_rows, file.path(out_dir, "cap-requirements.csv"),
          row.names = FALSE)
