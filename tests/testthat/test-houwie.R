context("test-houwie.R")

## hOUwie evaluates a stochastic likelihood, so most of what can go wrong shows up as a
## plausible looking number rather than an error. These check the invariants that pin
## that number down: a cached likelihood must belong to the parameters it is returned
## for, the answer must not depend on where the trait sits on the number line, and every
## documented optimizer must actually run.

make.test.simmap <- function(){
    data(tworegime)
    map <- getMapFromNode(tree, trait[,2], tree$node.label, 0.5)
    getMapFromSubstHistory(list(map), tree)[[1]]
}

test_that("the likelihood cache is keyed on the whole parameter vector", {
    skip_on_cran()

    ## Two regimes differing only in which of them carries the large sigma are
    ## different models. A cache that matched on the sum of the parameters rather than
    ## on the parameters themselves treated them as the same point and handed back the
    ## likelihood of whichever was evaluated first.
    data(tworegime)
    simmap <- make.test.simmap()
    hOUwie.dat <- organizeHOUwieDat(trait, "none", TRUE)
    index.disc <- getDiscreteModel(hOUwie.dat$data.cor, "ER", 1, FALSE, TRUE)
    index.disc[index.disc == 0] <- NA
    index.cont <- getOUParamStructure("OUMV", 2, 1, FALSE)
    all.paths <- lapply(1:(Nnode(simmap) + Ntip(simmap)), function(x) getPathToRoot(simmap, x))
    edge_liks_list <- getEdgeLiks(simmap, hOUwie.dat$data.cor, 2, 1, max(branching.times(simmap)) + 1)

    n_p <- max(index.disc, na.rm = TRUE) + max(index.cont, na.rm = TRUE)
    empty.cache <- function(){
        as.data.table(data.frame(matrix(c(0, rep(1e5, n_p)), byrow = TRUE, ncol = n_p + 1, nrow = 50)))
    }
    lik <- function(p, cache){
        hOUwie.fixed.dev(p = log(p), simmaps = list(simmap), data = hOUwie.dat$data.ou,
                         rate.cat = 1, tip.fog = "none", index.disc = index.disc,
                         index.cont = index.cont, root.p = "yang",
                         edge_liks_list = edge_liks_list, all.paths = all.paths,
                         split.liks = FALSE, global_liks_mat = cache)
    }

    # c(transition rate, alpha, sigma 1, sigma 2, theta 1, theta 2)
    pars    <- c(1, 2, 0.5, 4, 1, 3)
    swapped <- c(1, 2, 4, 0.5, 1, 3)

    shared <- empty.cache()
    first <- lik(pars, shared)
    # the same point must come back off the cache unchanged
    expect_equal(lik(pars, shared), first)
    # a permutation of it must not
    expect_equal(lik(swapped, shared), lik(swapped, empty.cache()))
    expect_false(isTRUE(all.equal(lik(swapped, empty.cache()), first)))
})

test_that("hOUwie does not depend on where the trait sits on the number line", {
    skip_on_cran()

    ## A trait with negative values is fit against a trait mean shifted right by 50 so
    ## that the thetas stay positive for a search on the log scale. Shifting the data by
    ## a constant and the thetas with it is the same model, so both the likelihood and
    ## the reconstruction have to be unchanged.
    data(tworegime)
    shifted <- trait
    shifted[,3] <- trait[,3] + 10

    pars         <- c(0.1, 1, 1, 0.5, 1.5)
    pars.shifted <- c(0.1, 1, 1, 10.5, 11.5)

    ## the map set is sampled, so both fits have to start from the same seed for the
    ## comparison to be about the shift rather than about the draw
    set.seed(1)
    fit <- hOUwie(tree, trait, rate.cat = 1, discrete_model = "ER",
                  continuous_model = "OUM", nSim = 10, quiet = TRUE, p = pars)
    set.seed(1)
    fit.shifted <- hOUwie(tree, shifted, rate.cat = 1, discrete_model = "ER",
                          continuous_model = "OUM", nSim = 10, quiet = TRUE, p = pars.shifted)
    expect_equal(fit$loglik, fit.shifted$loglik)

    ## the reconstruction is what silently degrades: log() of a theta reported back on
    ## the original scale is NaN, and every node comes out flat
    nodes <- (Ntip(tree) + 1):(Ntip(tree) + 3)
    set.seed(1)
    recon <- silence(hOUwie.recon(fit, nodes = nodes))
    set.seed(1)
    recon.shifted <- silence(hOUwie.recon(fit.shifted, nodes = nodes))
    expect_equal(recon, recon.shifted)
    expect_false(any(is.na(recon)))
})

test_that("every documented optimizer runs", {
    skip_on_cran()

    data(tworegime)
    simmap <- make.test.simmap()
    for(optimizer in c("nlopt_ln", "sann")){
        opts <- if(optimizer == "sann"){
            list(max.call = 10, smooth = FALSE)
        }else{
            list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = "10", "ftol_rel" = .Machine$double.eps^0.25)
        }
        fit <- silence(hOUwie.fixed(list(simmap), trait, rate.cat = 1, discrete_model = "ER",
                                    continuous_model = "OUM", quiet = TRUE,
                                    optimizer = optimizer, opts = opts))
        expect_true(is.finite(fit$loglik), info = optimizer)
    }
})
