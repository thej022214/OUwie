context("test-houwie-walk.R")

## hOUwie shifts a negative trait right before optimizing and getHouwieObj shifts the
## thetas back out, so houwie_obj$p is on the original trait scale while the bounds the
## fit was made under are on the shifted one. dent_propose rejects a proposal when any
## parameter is outside the bounds, so a scalar lower bound of zero applied to a negative
## theta rejects every proposal that leaves that theta alone, and the walk collapses onto
## the few proposals that move every parameter at once.

walk.fit <- function(trait_offset, p, nSim = 5){
    set.seed(42)
    phy <- rcoal(20)
    phy$edge.length <- phy$edge.length / max(branching.times(phy))
    states <- sample(c(1, 2), 20, replace = TRUE)
    x <- rnorm(20, ifelse(states == 1, 1, 2), 0.2) + trait_offset
    dat <- data.frame(sp = phy$tip.label, reg = states, x = x)
    hOUwie(phy = phy, data = dat, rate.cat = 1, discrete_model = "ER",
           continuous_model = "OUM", nSim = nSim, p = p, quiet = TRUE)
}

## the walk's objective, called exactly as hOUwie.walk calls it
walk_objective <- function(houwie_obj, par){
    withHOUwieSeed(houwie_obj$crn_seed,
                   hOUwie(p = par, phy = houwie_obj$phy, data = houwie_obj$data,
                          rate.cat = houwie_obj$rate.cat, tip.fog = houwie_obj$tip.fog,
                          discrete_model = houwie_obj$index.disc,
                          continuous_model = houwie_obj$index.cont,
                          root.p = houwie_obj$root.p, nSim = houwie_obj$nSim,
                          sample_tips = houwie_obj$sample_tips,
                          sample_nodes = houwie_obj$sample_nodes,
                          adaptive_sampling = houwie_obj$adaptive_sampling,
                          common_random_numbers = TRUE, quiet = TRUE)$loglik)
}

test_that("the walk evaluates the reported optimum without shifting the thetas twice", {
    skip_on_cran()

    fit <- walk.fit(-60, c(0.1, 48, 43, -59, -58))
    expect_true(fit$negative_values)
    expect_gt(fit$trait_shift, 0)

    at_best <- walk_objective(fit, fit$p)
    # the thetas hOUwie would be handed if they had already been shifted onto the scale
    # it shifts them onto itself
    shifted_par <- fit$p
    shifted_par[4:5] <- shifted_par[4:5] + fit$trait_shift
    at_shifted <- walk_objective(fit, shifted_par)

    # the objective resamples its maps, so it reproduces the fit only up to that sampling.
    # a second shift moves the optima a whole trait range away from the data
    expect_lt(abs(at_best - fit$loglik), 0.01 * abs(at_shifted - fit$loglik))
    # and hOUwie undoes its own shift exactly, so the reported parameters come back
    expect_equal(walk_objective(fit, fit$p), at_best)
})

test_that("a walk on a negative trait fit explores theta on the original trait scale", {
    skip_on_cran()

    fit <- walk.fit(-60, c(0.1, 48, 43, -59, -58))
    set.seed(7)
    walked <- suppressWarnings(hOUwie.walk(fit, nsteps = 20, print_freq = 1e6))

    thetas <- c("theta_1", "theta_2")
    widths <- walked$all_ranges["upper.CI", thetas] -
        walked$all_ranges["lower.CI", thetas]
    expect_true(all(widths > 0))
    # the reported interval is the scale the fit reports its thetas on
    expect_true(all(walked$all_ranges["best", thetas] < 0))

    # the constraint the fit was optimized under is theta >= -trait_shift, which is what
    # a lower bound of zero on the shifted scale means for the reported thetas
    sampled <- walked$results[, thetas]
    expect_true(all(sampled >= -getObjTraitShift(fit)))
    expect_true(all(walked$results[, c("rate_1", "alpha_1", "sigma2_1")] >= 0))

    # a walk that could not move theta would leave every proposal with all parameters at
    # a single recycled value
    all_equal_rows <- apply(walked$results[, -1], 1,
                            function(x) length(unique(x)) == 1)
    expect_lt(mean(all_equal_rows), 0.5)
})

test_that("a fit with no trait shift is walked under exactly the bounds it was given", {
    skip_on_cran()

    fit <- walk.fit(2, c(0.1, 48, 43, 3, 4))
    expect_false(isTRUE(fit$negative_values))
    expect_equal(getObjTraitShift(fit), 0)

    # with no shift the per parameter bounds are the scalar bound repeated, and
    # dent_propose compares bounds elementwise, so the proposals are the same draws
    old_params <- c(rate_1 = 0.1, alpha_1 = 48, sigma2_1 = 43, theta_1 = 3, theta_2 = 4)
    sd_vector <- 0.1 * abs(old_params)
    set.seed(21)
    scalar_bound <- t(replicate(100, dent_propose(old_params, lower_bound = 0,
                                                  upper_bound = Inf, sd = sd_vector)))
    set.seed(21)
    vector_bound <- t(replicate(100, dent_propose(old_params, lower_bound = rep(0, 5),
                                                  upper_bound = rep(Inf, 5),
                                                  sd = sd_vector)))
    expect_identical(scalar_bound, vector_bound)

    set.seed(7)
    walked <- suppressWarnings(hOUwie.walk(fit, nsteps = 20, print_freq = 1e6))
    expect_true(all(walked$results[, -1] >= 0))
    expect_true(all(is.finite(as.matrix(walked$results))))
})
