context("test-likelihood.R")

test_that("testing OU1 likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OU1", scaleHeight=TRUE, root.station=TRUE, shift.point=0.5, quiet=TRUE)
    comparison <- identical(round(ouwiefit$loglik,5), -22.54063)
    expect_true(comparison)
})


test_that("testing OUM likelihood", {
    skip_on_cran()

    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUM", scaleHeight=TRUE, root.station=TRUE, shift.point=0.5, quiet=TRUE)
    comparison <- identical(round(ouwiefit$loglik,5), -19.75473)
    expect_true(comparison)
})


test_that("testing OUMV likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUMV", scaleHeight=TRUE, root.station=TRUE, shift.point=0.5, quiet=TRUE)
    comparison <- identical(round(ouwiefit$loglik,5), -14.92345)
    expect_true(comparison)
})







