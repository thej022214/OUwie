#!/usr/bin/env Rscript

# Validate the deterministic mixture pruning against an exhaustive enumeration
# of hOUwie's own midpoint history model on a tree small enough to enumerate.

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

candidates <- data.frame(
  candidate = c("Generating parameters", "fitmultiOU estimate from post",
                "hOUwie estimate from post"),
  theta_a = c(5, 5.0209, 4.4566072),
  theta_b = c(10, 10.1985, 88.5531307),
  alpha = c(0.5, 0.4919, 0.0123027),
  sigma_sq = c(0.2, 0.1839, 0.7764039),
  q = c(0.3, 0.347199, 0.1187636),
  stringsAsFactors = FALSE
)

# the same fixture the enumeration in diagnose-history-integration.R uses
set.seed(96001)
small_phy <- phytools::pbtree(n = 6, scale = 10)
truth_Q <- matrix(c(-0.3, 0.3, 0.3, -0.3), 2, 2, byrow = TRUE)
small_simulation <- hOUwie.sim(small_phy, Q = truth_Q, root.freqs = c(1, 0),
                               alpha = c(0.5, 0.5), sigma.sq = c(0.2, 0.2),
                               theta0 = 5, theta = c(5, 10))
small_data <- small_simulation$data
small_phy <- ape::reorder.phylo(small_phy, "pruningwise")
small_data <- OUwie:::matchTipsAndData(small_phy$tip.label, small_data)
houwie_data <- OUwie:::organizeHOUwieDat(small_data, "none", TRUE)
tree_plan <- OUwie:::getHOUwieTreePlan(
  small_phy,
  OUwie:::getEdgeLiks(small_phy, houwie_data$data.cor, 2, 1,
                      max(ape::branching.times(small_phy)) + 1)
)
tip_states <- as.integer(houwie_data$data.cor[, 2])
names(tip_states) <- houwie_data$data.cor[, 1]
tip_values <- houwie_data$data.ou[, 3]
names(tip_values) <- houwie_data$data.ou[, 1]

internal_nodes <- (Ntip(small_phy) + 1):(Ntip(small_phy) + Nnode(small_phy))
assignments <- expand.grid(rep(list(1:2), length(internal_nodes)))

exact_enumeration <- function(candidate){
  q <- candidate$q
  Q <- matrix(c(-q, q, q, -q), 2, 2, byrow = TRUE)
  edge_transitions <- lapply(small_phy$edge.length,
                             function(t) expm::expm(Q * t))
  joint <- numeric(nrow(assignments))
  for(i in seq_len(nrow(assignments))){
    states <- integer(Ntip(small_phy) + Nnode(small_phy))
    states[seq_len(Ntip(small_phy))] <- tip_states[small_phy$tip.label]
    states[internal_nodes] <- as.integer(assignments[i, ])
    discrete <- log(0.5) + sum(vapply(seq_len(nrow(small_phy$edge)), function(e){
      log(edge_transitions[[e]][states[small_phy$edge[e, 1]],
                                states[small_phy$edge[e, 2]]])
    }, numeric(1)))
    midpoint_map <- mapply(function(a, d, len){
      setNames(rep(len / 2, 2), c(states[a], states[d]))
    }, a = small_phy$edge[, 1], d = small_phy$edge[, 2],
    len = small_phy$edge.length, SIMPLIFY = FALSE)
    continuous <- OUwie:::OUwie.basic(
      small_phy, houwie_data$data.ou, simmap.tree = TRUE, scaleHeight = FALSE,
      alpha = rep(candidate$alpha, 2), sigma.sq = rep(candidate$sigma_sq, 2),
      theta = unname(c(candidate$theta_a, candidate$theta_b)),
      algorithm = "three.point", tip.paths = tree_plan$tip.paths,
      tip.fog = "none", map = midpoint_map, tree.plan = tree_plan
    )
    joint[i] <- discrete + continuous
  }
  log_sum_exp(joint)
}

rows <- list()
for(i in seq_len(nrow(candidates))){
  candidate <- candidates[i, ]
  Q <- matrix(c(-candidate$q, candidate$q, candidate$q, -candidate$q), 2, 2,
              byrow = TRUE)
  alpha <- rep(candidate$alpha, 2)
  sigma_sq <- rep(candidate$sigma_sq, 2)
  theta <- c(candidate$theta_a, candidate$theta_b)
  exact <- exact_enumeration(candidate)
  for(cap in c(Inf, 32, 8, 4, 2, 1)){
    value <- mixture_pruning_loglik(
      small_phy, tip_states, tip_values, Q, alpha, sigma_sq, theta,
      root_p = c(0.5, 0.5), resolution = 1L, max_components = cap
    )
    rows[[length(rows) + 1L]] <- data.frame(
      candidate = candidate$candidate, max_components = cap,
      enumerated = exact, pruned = value, difference = value - exact
    )
  }
}
result <- do.call(rbind, rows)
print(result, row.names = FALSE, digits = 10)
