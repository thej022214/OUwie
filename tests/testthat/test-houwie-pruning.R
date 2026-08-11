context("test-houwie-pruning.R")

# Build a tree small enough that every assignment of regimes to internal nodes
# can be enumerated, so the pruned likelihood has an exact reference rather than
# another approximation to compare against.
pruning_fixture <- function(seed = 96001, n_tip = 6){
  set.seed(seed)
  phy <- phytools::pbtree(n = n_tip, scale = 10)
  Q <- matrix(c(-0.3, 0.3, 0.3, -0.3), 2, 2, byrow = TRUE)
  simulation <- hOUwie.sim(phy, Q = Q, root.freqs = c(1, 0),
                           alpha = c(0.5, 0.5), sigma.sq = c(0.2, 0.2),
                           theta0 = 5, theta = c(5, 10))
  phy <- ape::reorder.phylo(phy, "pruningwise")
  data <- OUwie:::matchTipsAndData(phy$tip.label, simulation$data)
  organized <- OUwie:::organizeHOUwieDat(data, "none", TRUE)
  edge_liks <- OUwie:::getEdgeLiks(phy, organized$data.cor, 2, 1,
                                   max(ape::branching.times(phy)) + 1)
  list(phy = phy, data = data, organized = organized, edge_liks = edge_liks,
       tree_plan = OUwie:::getHOUwieTreePlan(phy, edge_liks))
}

# The likelihood of hOUwie's own history model, summed over every assignment.
enumerate_midpoint_lik <- function(fixture, q, alpha, sigma.sq, theta){
  phy <- fixture$phy
  n_tip <- ape::Ntip(phy)
  tip_states <- as.integer(fixture$organized$data.cor[, 2])
  names(tip_states) <- fixture$organized$data.cor[, 1]
  tip_states <- tip_states[phy$tip.label]
  internal <- (n_tip + 1):(n_tip + ape::Nnode(phy))
  assignments <- expand.grid(rep(list(1:2), length(internal)))
  Q <- matrix(c(-q, q, q, -q), 2, 2, byrow = TRUE)
  transitions <- lapply(phy$edge.length, function(t) expm::expm(Q * t))
  terms <- vapply(seq_len(nrow(assignments)), function(i){
    states <- integer(n_tip + ape::Nnode(phy))
    states[seq_len(n_tip)] <- tip_states
    states[internal] <- as.integer(assignments[i, ])
    discrete <- log(0.5) + sum(vapply(seq_len(nrow(phy$edge)), function(e){
      log(transitions[[e]][states[phy$edge[e, 1]], states[phy$edge[e, 2]]])
    }, numeric(1)))
    midpoint <- mapply(function(a, d, len){
      setNames(rep(len / 2, 2), c(states[a], states[d]))
    }, a = phy$edge[, 1], d = phy$edge[, 2], len = phy$edge.length,
    SIMPLIFY = FALSE)
    discrete + OUwie:::OUwie.basic(
      phy, fixture$organized$data.ou, simmap.tree = TRUE, scaleHeight = FALSE,
      alpha = alpha, sigma.sq = sigma.sq, theta = theta,
      algorithm = "three.point", tip.paths = fixture$tree_plan$tip.paths,
      tip.fog = "none", map = midpoint, tree.plan = fixture$tree_plan)
  }, numeric(1))
  OUwie:::logSumExpSafe(terms)
}

pruned_lik <- function(fixture, q, alpha, sigma.sq, theta,
                       max_components = Inf){
  phy <- fixture$phy
  tip_values <- fixture$organized$data.ou[
    match(phy$tip.label, fixture$organized$data.ou[, 1]), 3]
  Q <- matrix(c(-q, q, q, -q), 2, 2, byrow = TRUE)
  as.numeric(OUwie:::houwiePruningLik(
    phy = phy,
    tip_state_sets = OUwie:::getTipStateSets(phy, fixture$edge_liks, 2),
    tip_values = tip_values, Q = Q, alpha = alpha, sigma.sq = sigma.sq,
    theta = theta, root.p = c(0.5, 0.5), resolution = 1L,
    max_components = max_components))
}

test_that("the pruned likelihood equals an exhaustive enumeration", {
  skip_on_cran()
  fixture <- pruning_fixture()
  settings <- list(
    list(q = 0.3, alpha = c(0.5, 0.5), sigma.sq = c(0.2, 0.2),
         theta = c(5, 10)),
    # widely separated optima, where merging components is most dangerous
    list(q = 0.1187636, alpha = c(0.5, 0.5), sigma.sq = c(0.7764039, 0.7764039),
         theta = c(4.4566072, 88.5531307)),
    # regime specific alpha and sigma, i.e. an OUMVA parameterization
    list(q = 0.25, alpha = c(0.30, 1.20), sigma.sq = c(0.15, 0.60),
         theta = c(4, 12)),
    # one regime at the Brownian limit
    list(q = 0.25, alpha = c(1e-10, 0.90), sigma.sq = c(0.25, 0.40),
         theta = c(5, 11))
  )
  for(setting in settings){
    exact <- enumerate_midpoint_lik(fixture, setting$q, setting$alpha,
                                    setting$sigma.sq, setting$theta)
    pruned <- pruned_lik(fixture, setting$q, setting$alpha, setting$sigma.sq,
                         setting$theta)
    expect_equal(pruned, exact, tolerance = 1e-6)
  }
})

test_that("the pruned likelihood is deterministic", {
  skip_on_cran()
  fixture <- pruning_fixture()
  first <- pruned_lik(fixture, 0.3, c(0.5, 0.5), c(0.2, 0.2), c(5, 10),
                      max_components = 4)
  second <- pruned_lik(fixture, 0.3, c(0.5, 0.5), c(0.2, 0.2), c(5, 10),
                       max_components = 4)
  expect_identical(first, second)
})

test_that("merging components never gains mass over the exact mixture", {
  skip_on_cran()
  fixture <- pruning_fixture()
  exact <- pruned_lik(fixture, 0.3, c(0.5, 0.5), c(0.2, 0.2), c(5, 10))
  for(cap in c(2, 4, 8)){
    reduced <- pruned_lik(fixture, 0.3, c(0.5, 0.5), c(0.2, 0.2), c(5, 10),
                          max_components = cap)
    expect_true(is.finite(reduced))
    expect_lt(abs(reduced - exact), 1)
  }
})

