context("test-anc-state.R")

test_that("ancestral state estimation works", {
  data("tworegime")
  fitted <- OUwie(tree,trait,model=c("OUMV"),root.station=TRUE)
  recon <- OUwie.anc(fitted)
  phy.stubbed <- attach.stub.taxa(fitted$phy)
  data.stubbed <- add.stub.taxa.to.data(phy.stubbed, fitted$data)

})
