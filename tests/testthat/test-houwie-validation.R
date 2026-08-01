context("test-houwie-validation.R")

## hOUwie accepts a lot of user input that it never checks, and most of the ways that
## input can be wrong end in a plausible looking number, a NULL return, or a subscript
## error deep inside a likelihood rather than a message about the input. These pin down
## the checks: the number of states has to be read as a number, a fit produced without
## sampled maps has to say so rather than reconstruct from nothing, a parameter vector of
## the wrong length has to be rejected whichever side of zero the trait sits on, and a
## tree tip with no data has to be named.

quick.opts <- function(maxeval = 2){
    list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = as.character(maxeval),
         "ftol_rel" = .Machine$double.eps^0.25)
}

make.many.state.data <- function(){
    set.seed(1)
    phy <- rcoal(24)
    phy$tip.label <- paste0("t", 1:24)
    states <- rep(1:12, 2)
    data <- data.frame(sp = phy$tip.label,
                       reg = as.character(states),
                       x = rnorm(24, states, 0.2))
    list(phy = phy, data = data)
}

test_that("the number of states is read from the state codes rather than their spelling", {
    expect_equal(getNumberOfStates(as.character(1:12)), 12)
    expect_equal(getNumberOfStates(c("1", "2&3", "10")), 10)
    # a dataset in which every tip is polymorphic leaves no plain codes at all
    expect_equal(getNumberOfStates(c("2&3", "1&10")), 10)
    expect_equal(getNumberOfStates(factor(c("9", "11"))), 11)
})

test_that("hOUwie.recon and hOUwie.fixed accept a dataset with more than nine states", {
    skip_on_cran()

    sim <- make.many.state.data()
    set.seed(2)
    fit <- hOUwie(sim$phy, sim$data, rate.cat = 1, discrete_model = "ER",
                  continuous_model = "OUM", nSim = 2, opts = quick.opts(),
                  quiet = TRUE)
    expect_equal(ncol(fit$index.cont), 12)

    recon <- hOUwie.recon(fit, nodes = length(sim$phy$tip.label) + 1)
    expect_equal(ncol(recon), 12)
    expect_equal(sum(recon), 1, tolerance = 1e-6)

    fixed_fit <- hOUwie.fixed(fit$simmaps[1], sim$data, rate.cat = 1,
                              discrete_model = "ER", continuous_model = "OUM",
                              opts = quick.opts(), quiet = TRUE, make_numeric = FALSE)
    expect_s3_class(fixed_fit, "houwie")
    expect_equal(ncol(fixed_fit$index.cont), 12)
    expect_equal(ncol(fixed_fit$solution.cont), 12)
})

test_that("every tip gets a starting regime when there are more than nine states", {
    skip_on_cran()

    sim <- make.many.state.data()
    # the state by rate category grid that assigns each tip a starting regime is built
    # from the number of states, and a tip left out of it keeps the 0 it was initialised
    # with, which getIP.theta then averages into no theta at all
    seen_states <- NULL
    reference <- getIP.theta
    local_mocked_bindings(
        getIP.theta = function(x, states, index){
            seen_states <<- states
            reference(x, states, index)
        }, .package = "OUwie")

    set.seed(7)
    fit <- hOUwie(sim$phy, sim$data, rate.cat = 2, discrete_model = "ER",
                  continuous_model = "OUM", nSim = 2, opts = quick.opts(),
                  quiet = TRUE)

    expect_s3_class(fit, "houwie")
    expect_equal(length(seen_states), length(sim$phy$tip.label))
    expect_true(all(seen_states > 0))
    expect_true(all(seen_states <= 12 * 2))
})

test_that("reconstruction refuses a fit that has no sampled maps", {
    skip_on_cran()

    data(tworegime)
    map <- getMapFromNode(tree, trait[,2], tree$node.label, 0.5)
    simmap <- getMapFromSubstHistory(list(map), tree)[[1]]
    fixed_fit <- hOUwie.fixed(list(simmap), trait, rate.cat = 1,
                              discrete_model = "ER", continuous_model = "OUM",
                              opts = quick.opts(), quiet = TRUE)

    expect_null(fixed_fit$nSim)
    expect_error(hOUwie.recon(fixed_fit, nodes = length(tree$tip.label) + 1),
                 "fit by hOUwie.fixed", fixed = TRUE)
})

test_that("hOUwie.fixed keeps its optimization when reconstruction is impossible", {
    skip_on_cran()

    data(tworegime)
    map <- getMapFromNode(tree, trait[,2], tree$node.label, 0.5)
    simmap <- getMapFromSubstHistory(list(map), tree)[[1]]
    fit_fixed <- function(){
        hOUwie.fixed(list(simmap), trait, rate.cat = 1, discrete_model = "ER",
                     continuous_model = "OUM", opts = quick.opts(), quiet = TRUE,
                     recon = TRUE, nodes = length(tree$tip.label) + 1)
    }

    # asking for a reconstruction is a request the fit cannot honour, so it has to
    # arrive as one warning rather than an error dumped to the console by try()
    warnings <- character(0)
    err <- textConnection("stderr_capture", "w", local = TRUE)
    sink(err, type = "message")
    fixed_fit <- withCallingHandlers(fit_fixed(),
        warning = function(w){
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    sink(type = "message")
    close(err)

    expect_s3_class(fixed_fit, "houwie")
    expect_true(is.finite(fixed_fit$loglik))
    expect_equal(length(fixed_fit$p), 5)
    expect_null(fixed_fit$recon)
    expect_equal(length(warnings), 1)
    expect_match(warnings, "hOUwie.fixed", fixed = TRUE)
    expect_false(any(grepl("Error", stderr_capture)))
})

test_that("a wrong length p is rejected whichever side of zero the trait sits on", {
    skip_on_cran()

    data(tworegime)
    positive <- negative <- trait
    positive[,3] <- positive[,3] - min(positive[,3]) + 1
    negative[,3] <- negative[,3] - max(negative[,3]) - 1

    fit_with <- function(dat, ...){
        hOUwie(tree, dat, rate.cat = 1, discrete_model = "ER",
               continuous_model = "OUM", nSim = 2, quiet = TRUE, ...)
    }

    expect_error(fit_with(positive, p = c(0.1, 0.1, 0.1)),
                 "You have supplied 3 in p, but the model structure requires 5",
                 fixed = TRUE)
    expect_error(fit_with(negative, p = c(0.1, 0.1, 0.1)),
                 "You have supplied 3 in p, but the model structure requires 5",
                 fixed = TRUE)
    expect_error(fit_with(negative, ip = c(0.1, 0.1, 0.1), opts = quick.opts()),
                 "You have supplied 3 in ip, but the model structure requires 5",
                 fixed = TRUE)
})

test_that("model averaging works when a model retained a single map", {
    skip_on_cran()

    data(tworegime)
    set.seed(3)
    one_map <- lapply(c("OUM", "BM1"), function(x)
        hOUwie(tree, trait, rate.cat = 1, discrete_model = "ER",
               continuous_model = x, nSim = 1, opts = quick.opts(), quiet = TRUE))
    expect_equal(unlist(lapply(one_map, function(x) length(x$simmaps))), c(1, 1))

    averaged <- getModelAvgParams(one_map)
    expect_equal(nrow(averaged), length(tree$tip.label))
    expect_true(all(is.finite(averaged$expected_mean)))
    expect_true(all(is.finite(averaged$expected_var)))
})

test_that("expected values over many maps still use the high mass maps", {
    skip_on_cran()

    data(tworegime)
    set.seed(4)
    many_maps <- hOUwie(tree, trait, rate.cat = 1, discrete_model = "ER",
                        continuous_model = "OUM", nSim = 10,
                        opts = quick.opts(), quiet = TRUE)
    expect_gt(length(many_maps$simmaps), 1)

    joint <- many_maps$all_cont_liks + many_maps$all_disc_liks
    weights <- exp(joint - max(joint))/sum(exp(joint - max(joint)))
    kept <- weights > 0.0099
    kept[2] <- TRUE
    weights <- weights[kept]/sum(weights[kept])
    expectations <- lapply(many_maps$simmaps[kept], function(x)
        getOUExpectations(x, many_maps$solution.cont))
    reference <- colSums(t(do.call(cbind, lapply(expectations, "[[", "expected_means"))) * weights)

    expected <- getExpectedValues(many_maps)
    expect_equal(expected$expected_mean,
                 unname(reference[1:length(tree$tip.label)]))
})

test_that("the model averaging convergence warning matches the gap it tests", {
    skip_on_cran()

    data(tworegime)
    set.seed(5)
    fits <- lapply(c("OUM", "BM1"), function(x)
        hOUwie(tree, trait, rate.cat = 1, discrete_model = "ER",
               continuous_model = x, nSim = 2, opts = quick.opts(), quiet = TRUE))
    fits[[2]]$BIC <- fits[[1]]$BIC + 2e5

    expect_warning(getModelAvgParams(fits), "exceeds 1e5", fixed = TRUE)
})

test_that("a tip with no data is named rather than silently dropped", {
    skip_on_cran()

    data(tworegime)
    missing_tip <- as.character(trait[1,1])
    incomplete <- trait[-1,]

    expect_error(hOUwie(tree, incomplete, rate.cat = 1, discrete_model = "ER",
                        continuous_model = "OUM", nSim = 2, opts = quick.opts(),
                        quiet = TRUE),
                 missing_tip, fixed = TRUE)

    map <- getMapFromNode(tree, trait[,2], tree$node.label, 0.5)
    simmap <- getMapFromSubstHistory(list(map), tree)[[1]]
    expect_error(hOUwie.fixed(list(simmap), incomplete, rate.cat = 1,
                              discrete_model = "ER", continuous_model = "OUM",
                              opts = quick.opts(), quiet = TRUE),
                 missing_tip, fixed = TRUE)
})

test_that("a data row with no tip is reported and removed", {
    skip_on_cran()

    data(tworegime)
    extra <- rbind(trait, trait[1,])
    extra[,1] <- as.character(extra[,1])
    extra[nrow(extra), 1] <- "not_a_tip"

    set.seed(6)
    fit <- expect_warning(hOUwie(tree, extra, rate.cat = 1, discrete_model = "ER",
                                 continuous_model = "OUM", nSim = 2,
                                 opts = quick.opts(), quiet = TRUE),
                          "not_a_tip", fixed = TRUE)
    expect_s3_class(fit, "houwie")
    expect_equal(nrow(fit$data), length(tree$tip.label))
})
