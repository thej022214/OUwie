context("test-houwie-shift.R")

## hOUwie optimizes on log(p), so every parameter it touches - the OU optima included -
## has to be strictly positive. A trait with negative values is therefore moved right
## before the fit and moved back before anything is reported. Nothing in that machinery
## raises an error when it goes wrong: log() of a negative number is a warning, and a
## trait moved by one amount against thetas moved by another is simply a different model.
## The invariants below are what a correct shift has to satisfy - it has to be large
## enough for the bounds it feeds, it has to leave the likelihood where it was, and it
## has to be undone in every place a fit reports the trait or its optima.

make.shift.tree <- function(){
    set.seed(31)
    phy <- rcoal(10)
    phy$edge.length <- phy$edge.length / max(branching.times(phy))
    phy
}

make.shift.data <- function(phy, centre){
    set.seed(19)
    data.frame(sp = phy$tip.label,
               reg = rep(c("1", "2"), length(phy$tip.label) / 2),
               x = rnorm(length(phy$tip.label), centre, 2))
}

test_that("the shift clears zero by a margin whatever the trait's magnitude", {
    skip_on_cran()

    ## the bounds on theta are min(trait)/10 and max(trait)*10, so a shifted minimum of
    ## exactly zero is no more usable than a negative one
    set.seed(11)
    for(centre in c(-1, -50, -200, -1e4)){
        x <- rnorm(20, centre, 2)
        expect_gt(min(x + getTraitShift(x)), 0)
    }

    ## an invariant negative trait has no spread to scale the margin by
    expect_gt(min(rep(-7, 5) + getTraitShift(rep(-7, 5))), 0)

    ## an all positive trait is left exactly where it is
    expect_identical(getTraitShift(c(0.5, 2, 9)), 0)
    expect_identical(getTraitShift(c(0, 3)), 0)
})

test_that("a trait far below zero can be fit at all", {
    skip_on_cran()

    phy <- make.shift.tree()
    dat <- make.shift.data(phy, -200)
    expect_lt(min(dat$x), -50)

    set.seed(55)
    fit <- hOUwie(phy, dat, rate.cat = 1, discrete_model = "ER",
                  continuous_model = "OU1", nSim = 5, quiet = TRUE)

    expect_true(is.finite(fit$loglik))
    expect_true(all(is.finite(fit$solution.cont[3,])))
    ## the optima are reported where the user's trait lives, not where it was fit
    expect_lt(max(fit$solution.cont[3,]), 0)
})

test_that("the likelihood does not depend on where the trait sits", {
    skip_on_cran()

    phy <- make.shift.tree()
    offset <- 250
    negative <- make.shift.data(phy, -200)
    positive <- negative
    positive$x <- negative$x + offset
    expect_lt(min(negative$x), -50)
    expect_gt(min(positive$x), 0)

    p <- c(0.5, 1.2, 0.9, -201, -198)
    set.seed(55)
    fit.negative <- hOUwie(phy, negative, rate.cat = 1, discrete_model = "ER",
                           continuous_model = "OUM", nSim = 10, p = p, quiet = TRUE)
    set.seed(55)
    fit.positive <- hOUwie(phy, positive, rate.cat = 1, discrete_model = "ER",
                           continuous_model = "OUM", nSim = 10,
                           p = c(p[1:3], p[4:5] + offset), quiet = TRUE)

    expect_true(fit.negative$negative_values)
    expect_false(fit.positive$negative_values)
    expect_equal(fit.negative$loglik, fit.positive$loglik)
    expect_equal(unname(fit.negative$solution.cont[3,]) + offset,
                 unname(fit.positive$solution.cont[3,]))

    ## a shift left half undone shows up here as a flat reconstruction on one side
    set.seed(77)
    recon.negative <- hOUwie.recon(fit.negative, nodes = "internal")
    set.seed(77)
    recon.positive <- hOUwie.recon(fit.positive, nodes = "internal")

    expect_equal(recon.negative, recon.positive)
    expect_gt(max(apply(recon.negative, 1, function(x) diff(range(x)))), 0.1)
})

test_that("a shifted fit reports the trait and its optima on one scale", {
    skip_on_cran()

    phy <- make.shift.tree()
    dat <- make.shift.data(phy, -200)

    set.seed(55)
    fit <- hOUwie(phy, dat, rate.cat = 1, discrete_model = "ER",
                  continuous_model = "OUM", nSim = 10,
                  p = c(0.5, 1.2, 0.9, -201, -198), quiet = TRUE)

    ## the reported trait, the copy the likelihood is recomputed against, and the optima
    ## all have to come back together - this is the pairing hOUwie.recon depends on
    expect_equal(fit$data[,3], dat$x)
    expect_equal(fit$hOUwie.dat$data.ou[,3], dat$x)
    expect_equal(unname(tail(fit$p, 2)), unname(fit$solution.cont[3,]))
    expect_equal(unname(fit$solution.cont[3,]), c(-201, -198))
})

test_that("a positive trait is fit exactly as if the shift did not exist", {
    skip_on_cran()

    phy <- make.shift.tree()
    dat <- make.shift.data(phy, 20)
    expect_gt(min(dat$x), 0)

    set.seed(101)
    fit <- hOUwie(phy, dat, rate.cat = 1, discrete_model = "ER",
                  continuous_model = "OUM", nSim = 10, quiet = TRUE)

    expect_false(fit$negative_values)
    expect_identical(fit$trait_shift, 0)
    expect_equal(fit$data[,3], dat$x)
    expect_equal(fit$hOUwie.dat$data.ou[,3], dat$x)
})

test_that("the shift finds the optima of a model without an alpha", {
    skip_on_cran()

    phy <- make.shift.tree()
    offset <- 250
    negative <- make.shift.data(phy, -200)
    positive <- negative
    positive$x <- negative$x + offset

    ## BM1 leaves the whole alpha row of index.cont as NA, so the optima sit at a
    ## different offset in p than they do under any OU model
    set.seed(9)
    fit.negative <- hOUwie(phy, negative, rate.cat = 1, discrete_model = "ER",
                           continuous_model = "BM1", nSim = 10,
                           p = c(0.5, 0.9, -199), quiet = TRUE)
    set.seed(9)
    fit.positive <- hOUwie(phy, positive, rate.cat = 1, discrete_model = "ER",
                           continuous_model = "BM1", nSim = 10,
                           p = c(0.5, 0.9, -199 + offset), quiet = TRUE)

    expect_true(all(is.na(fit.negative$index.cont[1,])))
    ## sigma sits next to the optima in p and must not have been moved with them
    expect_equal(unname(fit.negative$solution.cont[2,]), c(0.9, 0.9))
    expect_equal(unname(fit.negative$solution.cont[3,]), c(-199, -199))
    expect_equal(fit.negative$loglik, fit.positive$loglik)
    expect_equal(fit.negative$data[,3], negative$x)
})

test_that("only the trait mean moves when tip fog is supplied", {
    skip_on_cran()

    phy <- make.shift.tree()
    dat <- make.shift.data(phy, -200)
    dat$se <- rep(0.1, nrow(dat))

    set.seed(5)
    fit <- hOUwie(phy, dat, rate.cat = 1, discrete_model = "ER",
                  continuous_model = "OUM", nSim = 10, tip.fog = "known",
                  p = c(0.5, 1.2, 0.9, -201, -198), quiet = TRUE)

    expect_equal(fit$data[,3], dat$x)
    expect_equal(fit$data[,4], dat$se)
    expect_equal(fit$hOUwie.dat$data.ou[,3], dat$x)
    expect_equal(unname(fit$solution.cont[3,]), c(-201, -198))
})

test_that("a fit saved before the shift was recorded still reconstructs", {
    skip_on_cran()

    phy <- make.shift.tree()
    dat <- make.shift.data(phy, -5)

    set.seed(88)
    fit <- hOUwie(phy, dat, rate.cat = 1, discrete_model = "ER",
                  continuous_model = "OUM", nSim = 10,
                  p = c(0.6, 1.1, 0.8, -6, -4), quiet = TRUE)

    ## older fits carry the logical flag and no magnitude; 50 is the constant they were
    ## fit under and the only value that puts them back on the scale they were fit on
    stored <- fit
    stored$trait_shift <- NULL
    expect_equal(getObjTraitShift(stored), 50)
    expect_equal(getObjTraitShift(fit), fit$trait_shift)

    set.seed(9)
    recon.stored <- hOUwie.recon(stored, nodes = "internal")
    set.seed(9)
    recon.current <- hOUwie.recon(fit, nodes = "internal")

    expect_true(all(is.finite(recon.stored)))
    expect_equal(recon.stored, recon.current)
})
