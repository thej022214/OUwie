context("test-houwie-tree-plan.R")

make.plan.simmap <- function(){
    set.seed(73)
    phy <- rcoal(12)
    phy$edge.length <- phy$edge.length / max(branching.times(phy))
    tip.states <- rep(c(1, 3), length.out = Ntip(phy))
    node.states <- sample(c(1, 3), Nnode(phy), replace = TRUE)
    map <- getMapFromNode(phy, tip.states, node.states, 0.37)
    getMapFromSubstHistory(list(map), phy)[[1]]
}

reorder.plan.simmap <- function(simmap, index){
    simmap$edge <- simmap$edge[index, , drop = FALSE]
    simmap$edge.length <- simmap$edge.length[index]
    simmap$maps <- simmap$maps[index]
    simmap$mapped.edge <- simmap$mapped.edge[index, , drop = FALSE]
    attr(simmap, "order") <- NULL
    simmap
}

test_that("fused continuous moments match the established calculations", {
    simmap <- make.plan.simmap()
    set.seed(91)
    simmap <- reorder.plan.simmap(simmap, sample(seq_len(nrow(simmap$edge))))

    state.names <- colnames(simmap$mapped.edge)
    alpha <- c(0.4, 5, 1.2)
    sigma.sq <- c(0.7, 9, 1.8)
    theta <- c(-0.2, 100, 2.4)
    rate.mat <- rbind(alpha, sigma.sq, theta)[, as.integer(state.names), drop = FALSE]
    pars <- matrix(c(theta, sigma.sq, alpha), nrow = length(theta),
                   dimnames = list(as.character(seq_along(theta)),
                                   c("opt", "sig", "alp")))
    pars <- pars[as.integer(state.names), , drop = FALSE]
    edges <- cbind(seq_len(nrow(simmap$edge)), simmap$edge, 0, 0)
    root.edges <- which(simmap$edge[, 1] == Ntip(simmap) + 1L)
    root.state <- match(names(simmap$maps[[root.edges[2]]])[1], state.names)
    tip.paths <- lapply(seq_len(Ntip(simmap)), function(tip)
        getPathToRoot(simmap, tip))

    expected.weights <- weight.mat(simmap, edges, rate.mat,
                                   root.state = root.state,
                                   simmap.tree = TRUE,
                                   assume.station = TRUE)
    expected.transform <- transformPhy(simmap, simmap$maps, pars, tip.paths)
    actual <- continuousMapMoments(simmap, rate.mat, pars,
                                   root.state = root.state,
                                   assume.station = TRUE,
                                   map = simmap$maps,
                                   state.names = state.names)

    expect_equal(actual$W, expected.weights, tolerance = 1e-12)
    expect_equal(actual$tree$edge.length,
                 expected.transform$tree$edge.length, tolerance = 1e-12)
    expect_equal(actual$diag, expected.transform$diag, tolerance = 1e-12)
    expect_s3_class(actual$tree, "phylo")
    expect_null(actual$tree$maps)
    expect_null(actual$tree$mapped.edge)
})

test_that("compact and materialized map interfaces agree", {
    simmap <- reorder.phylo(make.plan.simmap(), "pruningwise")
    data <- data.frame(sp = simmap$tip.label,
                       reg = rep(c(1, 3), length.out = Ntip(simmap)),
                       x = seq(-0.5, 1.5, length.out = Ntip(simmap)))
    edge.liks <- getEdgeLiks(simmap, data[, 1:2], 3, 1,
                            max(branching.times(simmap)) + 1)
    plan <- getHOUwieTreePlan(simmap, edge.liks)
    args <- list(phy = simmap, data = data, simmap.tree = TRUE,
                 scaleHeight = FALSE, get.root.theta = TRUE,
                 alpha = c(0.4, 5, 1.2), sigma.sq = c(0.7, 9, 1.8),
                 theta = c(-0.2, 100, 2.4), algorithm = "three.point")

    materialized <- do.call(OUwie.basic, args)
    compact <- do.call(OUwie.basic,
                       c(args, list(map = simmap$maps,
                                    map.states = colnames(simmap$mapped.edge),
                                    tree.plan = plan)))
    materialized.expected <- do.call(OUwie.basic,
                                     c(args, list(return.expected.vals = TRUE)))
    compact.expected <- do.call(OUwie.basic,
                                c(args, list(map = simmap$maps,
                                             map.states = colnames(simmap$mapped.edge),
                                             tree.plan = plan,
                                             return.expected.vals = TRUE)))

    expect_equal(compact, materialized, tolerance = 1e-12)
    expect_equal(compact.expected, materialized.expected, tolerance = 1e-12)

    data$error <- seq(0.01, 0.12, length.out = Ntip(simmap))
    known.args <- args
    known.args$data <- data
    known.args$tip.fog <- "known"
    materialized.known <- do.call(OUwie.basic, known.args)
    compact.known <- do.call(OUwie.basic,
                             c(known.args, list(map = simmap$maps,
                                                map.states = colnames(simmap$mapped.edge),
                                                tree.plan = plan)))
    expect_equal(compact.known, materialized.known, tolerance = 1e-12)
})
