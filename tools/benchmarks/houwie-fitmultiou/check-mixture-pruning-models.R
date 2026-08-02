#!/usr/bin/env Rscript

# Does the mixture pruning cover the rest of hOUwie's model set? Two things are
# actually load-bearing and neither was exercised by the earlier checks:
#
#   1. alpha and sigma.sq differing between regimes (OUMVA, and every model it
#      subsumes - OUMA, OUMV, OUVA, OUA, OUV). Composing an edge then chains
#      intervals with different decay rates rather than a single one.
#   2. a tip compatible with more than one regime, which is what a hidden state
#      is: the observed state narrows the regime to a set, not to a value.
#
# Both are checked against exhaustive enumeration of the same midpoint history
# model, so the reference is exact rather than another approximation.

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

# enumerate every assignment of regimes to internal nodes, and to tips where the
# observation leaves a choice, under the one-switch-at-the-midpoint convention
enumerate_midpoint <- function(phy, tip_allowed, tip_values, Q, alpha,
                               sigma_sq, theta, root_p){
  n_tip <- Ntip(phy)
  n_node <- Nnode(phy)
  data_ou <- data.frame(sp = phy$tip.label, reg = 1,
                        x = unname(tip_values[phy$tip.label]),
                        stringsAsFactors = FALSE)
  tree_plan <- OUwie:::getHOUwieTreePlan(
    phy, OUwie:::getEdgeLiks(phy, data.frame(sp = phy$tip.label, reg = 1),
                             nrow(Q), 1, max(branching.times(phy)) + 1)
  )
  transitions <- lapply(phy$edge.length, function(t) expm::expm(Q * t))
  grid <- expand.grid(c(tip_allowed[phy$tip.label],
                        rep(list(seq_len(nrow(Q))), n_node)))
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
    continuous <- OUwie:::OUwie.basic(
      phy, data_ou, simmap.tree = TRUE, scaleHeight = FALSE,
      alpha = alpha, sigma.sq = sigma_sq, theta = theta,
      algorithm = "three.point", tip.paths = tree_plan$tip.paths,
      tip.fog = "none", map = midpoint, tree.plan = tree_plan,
      map.states = as.character(seq_len(nrow(Q)))
    )
    terms[i] <- discrete + continuous
  }
  log_sum_exp(terms)
}

report <- function(label, phy, tip_allowed, tip_values, Q, alpha, sigma_sq,
                   theta, root_p){
  exact <- enumerate_midpoint(phy, tip_allowed, tip_values, Q, alpha, sigma_sq,
                              theta, root_p)
  for(cap in c(Inf, 8, 4)){
    pruned <- mixture_pruning_loglik(phy, tip_allowed, tip_values, Q, alpha,
                                     sigma_sq, theta, root_p = root_p,
                                     resolution = 1L, max_components = cap)
    cat(sprintf("%-34s cap=%-4s enumerated=%13.7f  pruned=%13.7f  diff=%10.3e\n",
                label, ifelse(is.finite(cap), cap, "Inf"), exact, pruned,
                pruned - exact))
  }
}

# ---- 1. unequal alpha and sigma.sq between regimes (OUMVA) -----------------
set.seed(96001)
phy6 <- reorder.phylo(pbtree(n = 6, scale = 10), "pruningwise")
truth_Q <- matrix(c(-0.3, 0.3, 0.3, -0.3), 2, 2, byrow = TRUE)
sim <- hOUwie.sim(pbtree(n = 6, scale = 10), Q = truth_Q, root.freqs = c(1, 0),
                  alpha = c(0.5, 0.5), sigma.sq = c(0.2, 0.2), theta0 = 5,
                  theta = c(5, 10))
values6 <- setNames(sim$data[, 3], sim$data[, 1])[phy6$tip.label]
observed6 <- setNames(as.integer(sim$data[, 2]), sim$data[, 1])[phy6$tip.label]
q <- 0.25
Q2 <- matrix(c(-q, q, q, -q), 2, 2, byrow = TRUE)
allowed6 <- lapply(observed6, function(s) s)
names(allowed6) <- names(observed6)

report("OUMVA  unequal alpha+sigma+theta", phy6, allowed6, values6, Q2,
       alpha = c(0.30, 1.20), sigma_sq = c(0.15, 0.60),
       theta = c(4, 12), root_p = c(0.5, 0.5))
report("OUA    unequal alpha only", phy6, allowed6, values6, Q2,
       alpha = c(0.20, 1.50), sigma_sq = c(0.30, 0.30),
       theta = c(6, 6), root_p = c(0.5, 0.5))
report("BMS    alpha at the BM limit", phy6, allowed6, values6, Q2,
       alpha = c(1e-10, 1e-10), sigma_sq = c(0.20, 0.90),
       theta = c(5, 5), root_p = c(0.5, 0.5))
report("mixed  one BM regime one OU", phy6, allowed6, values6, Q2,
       alpha = c(1e-10, 0.90), sigma_sq = c(0.25, 0.40),
       theta = c(5, 11), root_p = c(0.5, 0.5))

# ---- 2. hidden states: four regimes, tips ambiguous within an observed state -
set.seed(4242)
phy4 <- reorder.phylo(pbtree(n = 4, scale = 6), "pruningwise")
values4 <- setNames(c(5.2, 6.1, 9.4, 10.2), phy4$tip.label)
# regimes 1,2 are the two rate categories of observed state a; 3,4 of state b
observed4 <- c(1L, 1L, 2L, 2L)
allowed4 <- lapply(observed4, function(s) if(s == 1L) c(1L, 2L) else c(3L, 4L))
names(allowed4) <- phy4$tip.label
rate_hidden <- 0.15
rate_observed <- 0.30
Q4 <- matrix(0, 4, 4)
Q4[1, 2] <- Q4[2, 1] <- Q4[3, 4] <- Q4[4, 3] <- rate_hidden
Q4[1, 3] <- Q4[3, 1] <- Q4[2, 4] <- Q4[4, 2] <- rate_observed
diag(Q4) <- -rowSums(Q4)

report("hidden 4 regimes, ambiguous tips", phy4, allowed4, values4, Q4,
       alpha = c(0.30, 0.90, 0.30, 0.90), sigma_sq = c(0.20, 0.50, 0.20, 0.50),
       theta = c(5, 5, 10, 10), root_p = rep(0.25, 4))
report("hidden 4 regimes, all differ", phy4, allowed4, values4, Q4,
       alpha = c(0.25, 1.10, 0.45, 0.80), sigma_sq = c(0.15, 0.55, 0.30, 0.70),
       theta = c(4.5, 6.0, 9.5, 11.5), root_p = rep(0.25, 4))
