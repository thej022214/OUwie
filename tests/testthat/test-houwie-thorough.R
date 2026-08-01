context("test-houwie-thorough.R")

## hOUwie.thorough restarts every model from a model averaged starting point, and
## runSingleThorough is the only place those averages are checked before they reach
## hOUwie.fixed. The search is over log(p), so a start has to be strictly positive on the
## scale it is finally logged on - and that scale is not the same for every parameter.
## hOUwie.fixed moves a negative trait right and carries the theta block of ip with it, so
## an optimum is only ever logged after the shift, while rates, alphas and sigmas are
## logged exactly as they are handed over. Judging thetas before the shift throws the
## averaging away for every negative trait; not judging them at all lets an optimum that
## is still non-positive after the shift reach log(), which is a warning and not an error.

thorough.cache <- new.env()

make.thorough.fit <- function(centre, tip.fog = "none"){
    key <- paste(centre, tip.fog, sep = "-")
    if(!is.null(thorough.cache[[key]])){
        return(thorough.cache[[key]])
    }
    set.seed(31)
    phy <- rcoal(8)
    phy$edge.length <- phy$edge.length / max(branching.times(phy))
    set.seed(19)
    dat <- data.frame(sp = phy$tip.label,
                      reg = rep(c("1", "2"), 4),
                      x = rnorm(8, centre, 2))
    if(tip.fog == "known"){
        dat$err <- rep(0.1, length(phy$tip.label))
    }
    set.seed(102)
    invisible(capture.output(
        fit <- hOUwie(phy, dat, rate.cat = 1, discrete_model = "ER",
                      continuous_model = "OUM", nSim = 3, quiet = TRUE,
                      tip.fog = tip.fog)))
    thorough.cache[[key]] <- fit
    fit
}

## the averaged tables runSingleThorough reduces to a starting vector. ER over two states
## gives one rate and OUM one alpha, one sigma and two optima, so ip is c(rate, alpha,
## sigma, theta, theta) and the optima are its last two entries
thorough.start <- function(fit, alpha, sigma, thetas, rate = 0.5){
    mod_avg_disc <- fit$solution.disc
    mod_avg_disc[!is.na(mod_avg_disc)] <- rate
    mod_avg_cont <- fit$solution.cont
    mod_avg_cont[1,] <- alpha
    mod_avg_cont[2,] <- sigma
    mod_avg_cont[3,] <- thetas
    list(mod_avg_disc = mod_avg_disc, mod_avg_cont = mod_avg_cont)
}

## the ip hOUwie.fixed records is the one it optimized from, so a NULL here is the
## fallback to default starting values and anything else is the averaged vector
run.thorough.start <- function(fit, init_pars){
    warned <- FALSE
    invisible(capture.output(withCallingHandlers(
        out <- runSingleThorough(fit, fit$simmaps[1:2], init_pars),
        warning = function(w){
            if(grepl("Model averaged", conditionMessage(w))){
                warned <<- TRUE
            }
            invokeRestart("muffleWarning")
        })))
    list(ip = out$ip, warned = warned)
}

test_that("optima averaged below zero are kept for a trait that lives below zero", {
    skip_on_cran()

    fit <- make.thorough.fit(-5)
    trait_shift <- getTraitShift(fit$hOUwie.dat$data.ou[,3])
    expect_gt(trait_shift, 0)

    res <- run.thorough.start(fit, thorough.start(fit, 0.5, 1, c(-5, -6)))

    expect_false(res$warned)
    expect_false(is.null(res$ip))
    ## the optima arrive shifted, everything else arrives as it was averaged
    expect_equal(res$ip, c(0.5, 0.5, 1, -5 + trait_shift, -6 + trait_shift))
})

test_that("a trait that lives above zero is averaged with no shift at all", {
    skip_on_cran()

    fit <- make.thorough.fit(10)
    expect_identical(getTraitShift(fit$hOUwie.dat$data.ou[,3]), 0)

    res <- run.thorough.start(fit, thorough.start(fit, 0.5, 1, c(9, 11)))

    expect_false(res$warned)
    expect_equal(res$ip, c(0.5, 0.5, 1, 9, 11))
})

test_that("optima that are still non-positive after the shift are rejected", {
    skip_on_cran()

    fit <- make.thorough.fit(-5)
    trait_shift <- getTraitShift(fit$hOUwie.dat$data.ou[,3])

    ## landing exactly on zero is no more usable than landing below it, log(0) is -Inf
    on_zero <- run.thorough.start(fit, thorough.start(fit, 0.5, 1, c(-trait_shift, -5)))
    expect_true(on_zero$warned)
    expect_null(on_zero$ip)

    below_zero <- run.thorough.start(fit, thorough.start(fit, 0.5, 1, c(-trait_shift - 1, -5)))
    expect_true(below_zero$warned)
    expect_null(below_zero$ip)

    ## an optimum that clears zero once shifted is usable however far below the data it is
    above_zero <- run.thorough.start(fit, thorough.start(fit, 0.5, 1, c(1e-3 - trait_shift, -5)))
    expect_false(above_zero$warned)
    expect_equal(above_zero$ip[4], 1e-3)
})

test_that("rates, alphas and sigmas are judged where they stand, not where thetas go", {
    skip_on_cran()

    ## the trait shift here is large enough to make any of these positive, and none of
    ## them is ever shifted, so all of them still have to be rejected
    fit <- make.thorough.fit(-5)
    expect_gt(getTraitShift(fit$hOUwie.dat$data.ou[,3]), 1)

    negative_alpha <- run.thorough.start(fit, thorough.start(fit, -1, 1, c(-5, -6)))
    expect_true(negative_alpha$warned)
    expect_null(negative_alpha$ip)

    negative_sigma <- run.thorough.start(fit, thorough.start(fit, 0.5, -1, c(-5, -6)))
    expect_true(negative_sigma$warned)
    expect_null(negative_sigma$ip)

    negative_rate <- run.thorough.start(fit, thorough.start(fit, 0.5, 1, c(-5, -6), rate = -1))
    expect_true(negative_rate$warned)
    expect_null(negative_rate$ip)

    zero_rate <- run.thorough.start(fit, thorough.start(fit, 0.5, 1, c(-5, -6), rate = 0))
    expect_true(zero_rate$warned)
    expect_null(zero_rate$ip)
})

test_that("averages that are not finite are rejected whichever process they belong to", {
    skip_on_cran()

    fit <- make.thorough.fit(-5)

    infinite_sigma <- run.thorough.start(fit, thorough.start(fit, 0.5, Inf, c(-5, -6)))
    expect_true(infinite_sigma$warned)
    expect_null(infinite_sigma$ip)

    undefined_theta <- run.thorough.start(fit, thorough.start(fit, 0.5, 1, c(NaN, -6)))
    expect_true(undefined_theta$warned)
    expect_null(undefined_theta$ip)

    missing_rate <- run.thorough.start(fit, thorough.start(fit, 0.5, 1, c(-5, -6), rate = NA))
    expect_true(missing_rate$warned)
    expect_null(missing_rate$ip)
})

test_that("the trait column is found the same way with tip fog as without", {
    skip_on_cran()

    ## with known tip fog the trait sits one column in from the end of data.ou, and the
    ## column past it is the fog. reading the last column would measure a shift of zero
    fit <- make.thorough.fit(-5, tip.fog = "known")
    expect_equal(dim(fit$hOUwie.dat$data.ou)[2], 4)
    trait_shift <- getTraitShift(fit$hOUwie.dat$data.ou[,3])
    expect_gt(trait_shift, 0)
    expect_identical(getTraitShift(fit$hOUwie.dat$data.ou[,4]), 0)

    res <- run.thorough.start(fit, thorough.start(fit, 0.5, 1, c(-5, -6)))

    expect_false(res$warned)
    expect_equal(res$ip, c(0.5, 0.5, 1, -5 + trait_shift, -6 + trait_shift))
})
