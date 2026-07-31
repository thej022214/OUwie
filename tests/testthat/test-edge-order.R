context("test-edge-order.R")

## The likelihood of a simmap is a property of the tree and its regime painting, not
## of the order the edges happen to be stored in. hOUwie hands OUwie.basic trees in
## pruningwise order while OUwie hands it cladewise ones, so both orders have to give
## the same answer.

reorder.simmap.for.test <- function(simmap, order){
    index <- reorder.phylo(simmap, order, index.only=TRUE)
    simmap$edge <- simmap$edge[index,]
    simmap$edge.length <- simmap$edge.length[index]
    simmap$maps <- simmap$maps[index]
    simmap$mapped.edge <- simmap$mapped.edge[index,]
    attr(simmap, "order") <- order
    simmap
}

make.test.simmap <- function(){
    data(tworegime)
    map <- getMapFromNode(tree, trait[,2], tree$node.label, 0.5)
    getMapFromSubstHistory(list(map), tree)[[1]]
}

test_that("OUwie.basic does not depend on edge ordering", {
    skip_on_cran()

    simmap <- make.test.simmap()
    data.ou <- data.frame(sp=trait[,1], reg=trait[,2], x=trait[,3])

    pars <- list(OU1   = list(alpha=c(2,2), sigma.sq=c(1,1),   theta=c(1.5,1.5)),
                 OUM   = list(alpha=c(2,2), sigma.sq=c(1,1),   theta=c(1,2)),
                 OUMV  = list(alpha=c(2,2), sigma.sq=c(1,0.5), theta=c(1,2)),
                 OUMA  = list(alpha=c(2,4), sigma.sq=c(1,1),   theta=c(1,2)),
                 OUMVA = list(alpha=c(2,4), sigma.sq=c(1,0.5), theta=c(1,2)))

    for(model in names(pars)){
        p <- pars[[model]]
        liks <- sapply(c("cladewise", "pruningwise"), function(order){
            phy <- reorder.simmap.for.test(simmap, order)
            tip.paths <- lapply(1:length(phy$tip.label), function(x) getPathToRoot(phy, x))
            OUwie.basic(phy, data.ou, simmap.tree=TRUE, alpha=p$alpha, sigma.sq=p$sigma.sq,
                        theta=p$theta, algorithm="three.point", tip.paths=tip.paths)
        })
        expect_equal(liks[["cladewise"]], liks[["pruningwise"]], info=model)
    }
})

test_that("OUwie.basic agrees with the explicit vcv likelihood", {
    skip_on_cran()

    ## OUwie.fixed's invert algorithm builds the OU vcv directly, so it is an
    ## independent check on the three.point route used by hOUwie.
    simmap <- make.test.simmap()
    data.ou <- data.frame(sp=trait[,1], reg=trait[,2], x=trait[,3])
    alpha <- c(2,4); sigma.sq <- c(1,0.5); theta <- c(1,2)

    phy <- reorder.simmap.for.test(simmap, "pruningwise")
    tip.paths <- lapply(1:length(phy$tip.label), function(x) getPathToRoot(phy, x))
    three.point <- OUwie.basic(phy, data.ou, simmap.tree=TRUE, alpha=alpha,
                               sigma.sq=sigma.sq, theta=theta, algorithm="three.point",
                               tip.paths=tip.paths)

    cladewise <- reorder.simmap.for.test(simmap, "cladewise")
    invert <- OUwie.fixed(cladewise, data.ou[match(cladewise$tip.label, data.ou$sp),],
                          model="OUMVA", simmap.tree=TRUE, algorithm="invert",
                          alpha=alpha, sigma.sq=sigma.sq, theta=theta,
                          check.identify=FALSE, quiet=TRUE)$loglik

    expect_equal(as.numeric(three.point), as.numeric(invert), tolerance=1e-6)
})
