context("test-houwie-guards.R")

## hOUwie.dev and hOUwie.fixed.dev return a deviance to the optimizer, and every way a
## likelihood evaluation can fail has to arrive there as the same finite 1e10 sentinel,
## because that is the only value hOUwie's multi-start filter recognises as a failure.
## These check that a genuine numerical failure produces it, and that a merely extreme
## but well defined likelihood does not.

houwie.dev.pieces <- function(continuous_model = "OUM"){
    data(tworegime)
    dat <- matchTipsAndData(tree$tip.label, trait)
    phy <- reorder.phylo(tree, "pruningwise")
    hdat <- organizeHOUwieDat(dat, "none", TRUE)
    nStates <- getNumberOfStates(hdat$data.cor[,2])
    index.disc <- getDiscreteModel(dat[,1:2], "ER", 1, FALSE, TRUE)
    index.disc[index.disc == 0] <- NA
    list(phy = phy,
         data = hdat$data.ou,
         index.disc = index.disc,
         index.cont = getOUParamStructure(continuous_model, nStates, 1, FALSE),
         edge_liks_list = getEdgeLiks(phy, hdat$data.cor, nStates, 1,
                                      max(branching.times(phy)) + 1),
         all.paths = lapply(1:(Nnode(phy) + Ntip(phy)),
                            function(x) getPathToRoot(phy, x)),
         simmaps = getMapFromSubstHistory(
             list(getMapFromNode(phy, dat[,2], phy$node.label, 0.5)), phy))
}

call.houwie.dev <- function(pieces, p, split.liks = FALSE, ...){
    hOUwie.dev(p = log(p), phy = pieces$phy, data = pieces$data, rate.cat = 1,
               tip.fog = "none", index.disc = pieces$index.disc,
               index.cont = pieces$index.cont, root.p = "yang",
               edge_liks_list = pieces$edge_liks_list, nSim = 5,
               all.paths = pieces$all.paths, split.liks = split.liks, ...)
}

call.houwie.fixed.dev <- function(pieces, p, split.liks = FALSE, root.p = "yang", ...){
    hOUwie.fixed.dev(p = log(p), simmaps = pieces$simmaps, data = pieces$data,
                     rate.cat = 1, tip.fog = "none", index.disc = pieces$index.disc,
                     index.cont = pieces$index.cont, root.p = root.p,
                     edge_liks_list = pieces$edge_liks_list,
                     all.paths = pieces$all.paths, split.liks = split.liks, ...)
}

## a tree whose tips are all but one in the same state, so that every map the sampler
## draws roots in that state. a root prior that gives that state no probability then
## leaves every map with probability zero.
one.state.pieces <- function(){
    set.seed(5)
    phy <- rcoal(12)
    phy$edge.length <- phy$edge.length / max(branching.times(phy))
    dat <- data.frame(sp = phy$tip.label, reg = c(rep(1, 11), 2),
                      x = rnorm(12, 1, 0.2))
    dat <- matchTipsAndData(phy$tip.label, dat)
    phy <- reorder.phylo(phy, "pruningwise")
    hdat <- organizeHOUwieDat(dat, "none", TRUE)
    nStates <- getNumberOfStates(hdat$data.cor[,2])
    index.disc <- getDiscreteModel(dat[,1:2], "ER", 1, FALSE, TRUE)
    index.disc[index.disc == 0] <- NA
    list(phy = phy,
         data = hdat$data.ou,
         index.disc = index.disc,
         index.cont = getOUParamStructure("OUM", nStates, 1, FALSE),
         edge_liks_list = getEdgeLiks(phy, hdat$data.cor, nStates, 1,
                                      max(branching.times(phy)) + 1),
         all.paths = lapply(1:(Nnode(phy) + Ntip(phy)),
                            function(x) getPathToRoot(phy, x)))
}

test_that("an underflowing set of maps returns the failure sentinel, not NaN", {
    skip_on_cran()

    pieces <- houwie.dev.pieces()
    # every map's continuous likelihood underflows to -Inf here, which leaves the
    # log sum exp over the maps as NaN
    underflow <- c(0.01, 1e-10, 1e-10, 1e200, 1e200)

    set.seed(3)
    dev <- call.houwie.dev(pieces, underflow, global_liks_mat = NULL)
    expect_true(is.finite(dev))
    expect_equal(dev, 1e10)

    set.seed(3)
    fixed_dev <- call.houwie.fixed.dev(pieces, underflow, global_liks_mat = NULL)
    expect_true(is.finite(fixed_dev))
    expect_equal(fixed_dev, 1e10)
})

test_that("the split likelihood report is a list even where the deviance is the sentinel", {
    skip_on_cran()

    pieces <- houwie.dev.pieces()
    underflow <- c(0.01, 1e-10, 1e-10, 1e200, 1e200)

    # hOUwie and hOUwie.fixed build their output object by indexing this report, so the
    # scalar sentinel can only ever be handed back to the optimizer
    set.seed(3)
    split <- call.houwie.dev(pieces, underflow, split.liks = TRUE,
                             global_liks_mat = NULL)
    expect_type(split, "list")
    expect_true(all(c("TotalLik", "DiscLik", "ContLik") %in% names(split)))

    set.seed(3)
    fixed_split <- call.houwie.fixed.dev(pieces, underflow, split.liks = TRUE,
                                         global_liks_mat = NULL)
    expect_type(fixed_split, "list")
    expect_true(all(c("TotalLik", "DiscLik", "ContLik") %in% names(fixed_split)))
})

test_that("the failure sentinel is what hOUwie treats as a failed optimization", {
    skip_on_cran()

    pieces <- houwie.dev.pieces()
    set.seed(3)
    dev <- call.houwie.dev(pieces, c(0.01, 1e-10, 1e-10, 1e200, 1e200),
                           global_liks_mat = NULL)
    expect_true(!is.finite(dev) | dev >= 1e10)
})

test_that("an extreme but well defined likelihood is reported rather than penalised", {
    skip_on_cran()

    pieces <- houwie.dev.pieces()
    # sigma.sq at its own lower bound is a legal point of the optimizer's parameter
    # box. the deviance there is enormous but finite, and the optimizer needs the
    # real value to be able to move away from it
    extreme <- c(0.01, 1e-10, 1e-10, 1.0, 2.0)

    set.seed(3)
    dev <- call.houwie.dev(pieces, extreme, global_liks_mat = NULL)
    expect_true(is.finite(dev))
    expect_true(dev > 1e10)

    set.seed(3)
    fixed_dev <- call.houwie.fixed.dev(pieces, extreme, global_liks_mat = NULL)
    expect_true(is.finite(fixed_dev))
    expect_true(fixed_dev > 1e10)
})

test_that("a draw with no usable map is a penalty to the optimizer and an error to a reporting caller", {
    skip_on_cran()

    pieces <- one.state.pieces()
    excluding_prior <- c(0, 1)
    p <- c(0.001, 1, 1, 1, 2)

    # hOUwie.dev draws its own maps, so an all-zero draw is transient and the optimizer
    # needs a finite value it can move away from
    set.seed(3)
    dev <- hOUwie.dev(p = log(p), phy = pieces$phy, data = pieces$data, rate.cat = 1,
                      tip.fog = "none", index.disc = pieces$index.disc,
                      index.cont = pieces$index.cont, root.p = excluding_prior,
                      edge_liks_list = pieces$edge_liks_list, nSim = 5,
                      all.paths = pieces$all.paths, split.liks = FALSE,
                      global_liks_mat = NULL)
    expect_equal(dev, 1e10)

    # the split likelihood report is handed to getHouwieObj, which indexes it as a list,
    # so there is nothing it can be given here
    set.seed(3)
    expect_error(
        hOUwie.dev(p = log(p), phy = pieces$phy, data = pieces$data, rate.cat = 1,
                   tip.fog = "none", index.disc = pieces$index.disc,
                   index.cont = pieces$index.cont, root.p = excluding_prior,
                   edge_liks_list = pieces$edge_liks_list, nSim = 5,
                   all.paths = pieces$all.paths, split.liks = TRUE,
                   global_liks_mat = NULL),
        "cannot be reported")
})

test_that("fixed maps that all have probability zero are reported as a modelling error", {
    skip_on_cran()

    pieces <- houwie.dev.pieces()
    # the supplied map has the root in state 1, and this prior gives that state no
    # probability. the maps are never resampled, so no parameter value recovers it
    excluding_prior <- c(0, 1)
    p <- c(0.01, 1, 1, 1, 2)

    expect_error(call.houwie.fixed.dev(pieces, p, split.liks = FALSE,
                                       root.p = excluding_prior,
                                       global_liks_mat = NULL),
                 "probability zero")
    expect_error(call.houwie.fixed.dev(pieces, p, split.liks = TRUE,
                                       root.p = excluding_prior,
                                       global_liks_mat = NULL),
                 "root.p")
})

test_that("a fixed fit under an excluding root prior names the problem instead of indexing a scalar", {
    skip_on_cran()

    data(tworegime)
    phy <- reorder.phylo(tree, "pruningwise")
    maps <- getMapFromSubstHistory(
        list(getMapFromNode(phy, trait[,2], tree$node.label, 0.5)), phy)

    expect_error(hOUwie.fixed(simmaps = maps, data = trait, rate.cat = 1,
                              discrete_model = "ER", continuous_model = "OUM",
                              root.p = c(0, 1), p = c(0.1, 1, 0.5, 2, 3),
                              quiet = TRUE),
                 "probability zero")
})

test_that("caching does not decide whether an evaluation is treated as a failure", {
    skip_on_cran()

    pieces <- houwie.dev.pieces()
    extreme <- c(0.01, 1e-10, 1e-10, 1.0, 2.0)
    liks_mat <- data.table::as.data.table(
        data.frame(matrix(c(0, rep(1e5, length(extreme))), byrow = TRUE,
                          ncol = length(extreme) + 1, nrow = 10)))

    set.seed(3)
    uncached <- call.houwie.dev(pieces, extreme, global_liks_mat = NULL)
    set.seed(3)
    cached <- call.houwie.dev(pieces, extreme, global_liks_mat = liks_mat)

    expect_equal(cached, uncached)
    expect_equal(as.numeric(liks_mat[1, 1]), -uncached)
})
