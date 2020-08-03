context("test-restart.R")

test_that("restarting works", {
    skip_on_cran()
  data("tworegime")
  fitted1 <- OUwie(tree,trait,model=c("OUMV"),root.station=FALSE, opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="5", "ftol_abs"=.1))
  fitted2 <- OUwie(tree,trait,model=c("OUMV"),root.station=FALSE, starting.vals=fitted1$new.start)
  testthat::expect_lt(fitted2$AIC, fitted1$AIC)
})
