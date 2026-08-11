context("test-model-averaging.R")

make.model.average.fit <- function(score, id = NULL){
    fit <- list(
        id = id,
        data = data.frame(sp = "sp1"),
        param.count = 1,
        phy = list(tip.label = "sp1"),
        AIC = score,
        AICc = score,
        BIC = score,
        loglik = -score / 2,
        DiscLik = -1,
        ContLik = -1,
        rate.cat = 1,
        index.disc = matrix(1, 1, 1),
        hOUwie.dat = list(
            PossibleTraits = "state1",
            data.cor = data.frame(sp = "sp1", state = 1)
        )
    )
    class(fit) <- "houwie"
    fit
}

mock.tip.values <- function(model){
    data.frame(
        rates = 1,
        alpha = 1,
        sigma.sq = 1,
        theta = 1,
        expected_mean = 1,
        expected_var = 1,
        row.names = "sp1"
    )
}

test_that("getModelTable diagnoses a one-model input directly", {
    expect_error(
        getModelTable(list(make.model.average.fit(0))),
        "Two or more models are needed"
    )
})

test_that("getModelAvgParams does not print intermediate results", {
    local_mocked_bindings(get_tip_values = mock.tip.values, .package = "OUwie")
    models <- list(make.model.average.fit(0), make.model.average.fit(10))

    expect_silent(result <- getModelAvgParams(models))
    expect_s3_class(result, "data.frame")
})

test_that("failed models are pruned relative to the best model", {
    used <- character()
    local_mocked_bindings(
        get_tip_values = function(model){
            used <<- c(used, model$id)
            mock.tip.values(model)
        },
        .package = "OUwie"
    )
    models <- list(
        best = make.model.average.fit(0, "best"),
        plausible = make.model.average.fit(100, "plausible"),
        failed = make.model.average.fit(200001, "failed")
    )

    expect_silent(getModelAvgParams(models, force = FALSE))
    expect_equal(used, c("best", "plausible"))
})

test_that("the convergence cutoff is configurable and validated", {
    used <- character()
    local_mocked_bindings(
        get_tip_values = function(model){
            used <<- c(used, model$id)
            mock.tip.values(model)
        },
        .package = "OUwie"
    )
    models <- list(
        best = make.model.average.fit(0, "best"),
        plausible = make.model.average.fit(25, "plausible"),
        failed = make.model.average.fit(100, "failed")
    )

    expect_silent(getModelAvgParams(models, force = FALSE, convergence_cutoff = 50))
    expect_equal(used, c("best", "plausible"))
    expect_error(
        getModelAvgParams(models, convergence_cutoff = 0),
        "single positive finite number"
    )
})

test_that("missing rows introduced while matching tips are reported early", {
    expect_error(
        OUwie:::matchTipsAndData(c("sp1", "sp2"), data.frame(sp = c("sp1", NA))),
        "Some tips in your phylogeny have no matching row in your data: sp2"
    )
})
