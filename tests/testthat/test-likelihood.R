context("test-likelihood.R")

test_that("testing OU1 likelihood stationary", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OU1", scaleHeight=TRUE, root.station=TRUE, shift.point=0.5, quiet=TRUE)
    comparison <- identical(round(ouwiefit$loglik,5), -22.54063)
    expect_true(comparison)
})


test_that("testing OUM likelihood stationary", {
    skip_on_cran()

    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUM", scaleHeight=TRUE, root.station=TRUE, shift.point=0.5, quiet=TRUE)
    comparison <- identical(round(ouwiefit$loglik,5), -19.75473)
    expect_true(comparison)
})

test_that("testing OU1 likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OU1", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, quiet=TRUE)
    comparison <- identical(round(ouwiefit$loglik,5), -21.74538)
    expect_true(comparison)
})


test_that("testing OUM likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUM", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, quiet=TRUE)
    comparison <- identical(round(ouwiefit$loglik,5), -19.51361)
    expect_true(comparison)
})

test_that("testing OUMV likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUMV", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, quiet=TRUE)
    comparison <- identical(round(ouwiefit$loglik,5), -14.79506)
    expect_true(comparison)
})

test_that("testing OUMA likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUMA", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, quiet=TRUE)
    comparison <- identical(round(ouwiefit$loglik,5), -19.42735)
    expect_true(comparison)
})


test_that("testing OUMVA likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUMVA", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, quiet=TRUE)
    comparison <- identical(round(ouwiefit$loglik,5), -13.97626)
    expect_true(comparison)
})


test_that("testing simmap", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(8)
    library(phytools)
    regs <- setNames(trait[,2], trait[,1])
    test <- make.simmap(tree, regs, model="ER")
    for(i in 1:dim(test$mapped.edge)[1]){
        entries <- test$mapped.edge[i,which(test$mapped.edge[i,] > 0)]
        test$mapped.edge[i,which(test$mapped.edge[i,] > 0)] <- sum(entries)/length(entries)
        maps <- test$maps[[i]]
        test$maps[[i]] <- rep(sum(maps)/length(maps), length(maps))
        names(test$maps[[i]]) <- names(maps)
    }
    ouwiefit.nodes <- OUwie(tree, trait, model="OUM", root.station=FALSE, shift.point=0.5, quiet=TRUE)
    ouwiefit.simmap <- OUwie(test, trait, model="OUM", simmap.tree=TRUE, root.station=FALSE, shift.point=0.5, quiet=TRUE)
    comparison <- identical(round(ouwiefit.nodes$loglik,5), round(ouwiefit.simmap$loglik,5))
    expect_true(comparison)
})


## For testing BM1 and BMS models:
#test_that("testing BM1", {
#    skip_on_cran()

#    library(phytools)
#    library(geiger)
#    library(OUwie)

    ## simulate some data
#    tree<-pbtree(n=26)
#    Q<-matrix(c(-1,1,1,-1),2,2,dimnames=list(letters[1:2],
#    letters[1:2]))
#    tree<-sim.history(tree,Q)
#    plot(tree)
#    x<-as.factor(getStates(tree,"tips"))
#    y<-fastBM(tree)

    ## fit using brownie.lite
#    brownie.lite(tree,y)

    ## fit using OUwie
#    test.data<-data.frame(Genus_species=tree$tip.label,Reg=x,X=y)
#    OUwie(tree,test.data,model="BM1",simmap.tree=TRUE, root.station=FALSE)

    ## fit using fitContinuous
#    fitContinuous(tree,y)
#})






