context("test-likelihood.R")

test_that("testing OU1 likelihood stationary", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OU1", scaleHeight=TRUE, root.station=TRUE, shift.point=0.5, algorithm="invert", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,5), -22.54063)
    expect_true(comparison)
})


test_that("testing OUM likelihood stationary", {
    skip_on_cran()

    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUM", scaleHeight=TRUE, root.station=TRUE, shift.point=0.5, algorithm="invert", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,5), -19.75473)
    expect_true(comparison)
})


test_that("testing BM1 likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="BM1", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="invert", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,5), -21.95911)
    expect_true(comparison)
})


test_that("testing BMS likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="BMS", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="invert", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,5), -17.85074)
    expect_true(comparison)
})


test_that("testing OU1 likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OU1", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="invert", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,5), -21.74538)
    expect_true(comparison)
})


test_that("testing OUM likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUM", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="invert", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,5), -19.51361)
    expect_true(comparison)
})


test_that("testing OUMV likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUMV", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="invert", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,5), -14.79506)
    expect_true(comparison)
})


test_that("testing OUMA likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUMA", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="invert", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,5), -19.42795)
    expect_true(comparison)
})


test_that("testing OUMVA likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUMVA", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="invert", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,5), -13.97626)
    expect_true(comparison)
})


#test_that("testing simmap", {
#    skip_on_cran()
    
#    data(tworegime)
#    set.seed(8)
#    library(phytools)
#    library(corHMM)
#    regs <- setNames(trait[,2], trait[,1])
#    test <- make.simmap(tree, regs, model="ER")
#    tree.corm <- reorder(tree, "pruningwise")
#    corm <- makeSimmap(tree.corm, trait[c(1,2)], test$Q, 1)[[1]]
#    for(i in 1:dim(test$mapped.edge)[1]){
#        entries <- test$mapped.edge[i,which(test$mapped.edge[i,] > 0)]
#        test$mapped.edge[i,which(test$mapped.edge[i,] > 0)] <- sum(entries)/length(entries)
#        maps <- test$maps[[i]]
#        test$maps[[i]] <- rep(sum(maps)/length(maps), length(maps))
#        names(test$maps[[i]]) <- names(maps)
#        entries <- corm$mapped.edge[i,which(corm$mapped.edge[i,] > 0)]
#        corm$mapped.edge[i,which(corm$mapped.edge[i,] > 0)] <- sum(entries)/length(entries)
#        maps <- corm$maps[[i]]
#        corm$maps[[i]] <- rep(sum(maps)/length(maps), length(maps))
#        names(corm$maps[[i]]) <- names(maps)
#    }
#    ouwiefit.nodes <- OUwie(tree, trait, model="OUM", root.station=FALSE, shift.point=0.5, algorithm="invert", quiet=TRUE)
#    ouwiefit.simmap <- OUwie(test, trait, model="OUM", simmap.tree=TRUE, root.station=FALSE, algorithm="invert", shift.point=0.5, quiet=TRUE)
#    ouwiefit.corm <- OUwie(corm, trait, model="OUM", simmap.tree=TRUE, root.station=FALSE, algorithm="invert", shift.point=0.5, quiet=TRUE)
#    # ouwiefit.3pta <- OUwie(test, trait, model="OUM", simmap.tree=TRUE, root.station=FALSE, algorithm="three.point", shift.point=0.5, quiet=TRUE)
#    # ouwiefit.3ptb <- OUwie(tree, trait, model="OUM", root.station=FALSE, shift.point=0.5, algorithm="three.point", quiet=TRUE)
#        comparison <- identical(round(ouwiefit.nodes$loglik,5), round(ouwiefit.simmap$loglik,5), round(ouwiefit.corm$loglik,5))
#    expect_true(comparison)
#})


test_that("testing BM1 likelihood invert vs three.point", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    
    sigma.sq=c(0.4669113, 0.4669113)
    theta=c(1.326483, 1.326483)
    
    BM1Invert <- OUwie.fixed(tree, trait, model=c("BM1"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, sigma.sq=sigma.sq, shift.point=0.5, algorithm="invert", check.identify=FALSE)
    BM13Point <- OUwie.fixed(tree, trait, model=c("BM1"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, sigma.sq=sigma.sq, theta=theta, shift.point=0.5, algorithm="three.point", check.identify=FALSE)
    comparison <- identical(round(as.numeric(BM1Invert$loglik),5), round(as.numeric(BM13Point$loglik),5))
    expect_true(comparison)
})


test_that("testing BMS likelihood invert vs three.point", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    
    sigma.sq=c(0.2424788, 0.7007112)
    theta <- c(1.4719377, 1.4719377)
    BMSInvert <- OUwie.fixed(tree, trait, model=c("BMS"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, sigma.sq=sigma.sq, shift.point=0.5, algorithm="invert", check.identify=FALSE)
    BMS3Point <- OUwie.fixed(tree, trait, model=c("BM1"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, sigma.sq=sigma.sq, theta=theta, shift.point=0.5, algorithm="three.point", check.identify=FALSE)
    comparison <- identical(round(as.numeric(BMSInvert$loglik),5), round(as.numeric(BMS3Point$loglik),5))
    expect_true(comparison)
})


test_that("testing OU1 likelihood invert vs three.point", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    
    alpha=c(0.358939, 0.3589399)
    sigma.sq=c(0.5197486, 0.5197486)
    theta=c( 1.3301447, 1.3301447)
    
    OU1Invert <- OUwie.fixed(tree, trait, model=c("OU1"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq,theta=theta, shift.point=0.5, algorithm="invert", check.identify=FALSE)
    OU13Point <- OUwie.fixed(tree, trait, model=c("OU1"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq,theta=theta, shift.point=0.5, algorithm="three.point", check.identify=FALSE)
    comparison <- identical(round(as.numeric(OU1Invert$loglik),5), round(as.numeric(OU13Point$loglik),5))
    expect_true(comparison)
})


test_that("testing OUM likelihood invert vs three.point", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    
    alpha=c(1.3916589, 1.3916589)
    sigma.sq=c(0.6545502, 0.6545502)
    theta=c(1.6751330, 0.4424138)

    OUMInvert <- OUwie.fixed(tree, trait, model=c("OUM"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq,theta=theta, shift.point=0.5, algorithm="invert", check.identify=FALSE)
    OUM3Point <- OUwie.fixed(tree, trait, model=c("OUM"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq,theta=theta, shift.point=0.5, algorithm="three.point", check.identify=FALSE)
    comparison <- identical(round(as.numeric(OUMInvert$loglik),5), round(as.numeric(OUM3Point$loglik),5))
    expect_true(comparison)
})


test_that("testing OUMV likelihood invert vs three.point", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    
    alpha=c(1.7110818, 1.711082)
    sigma.sq=c(0.3517019, 1.076479)
    theta=c(1.676894, 0.5563541)
    
    OUMVInvert <- OUwie.fixed(tree, trait, model=c("OUMV"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq,theta=theta, shift.point=0.5, algorithm="invert", check.identify=FALSE)
    OUMV3Point <- OUwie.fixed(tree, trait, model=c("OUMV"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq,theta=theta, shift.point=0.5, algorithm="three.point", check.identify=FALSE)
    comparison <- identical(round(as.numeric(OUMVInvert$loglik),5), round(as.numeric(OUMV3Point$loglik),5))
		
    expect_true(comparison)
})


test_that("testing OUMA likelihood invert vs three.point", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    
    alpha=c(1.6501816, 1.0294487)
    sigma.sq=c(0.7082462, 0.7082462)
    theta=c(1.6765718, 0.1516105)
    
    OUMAInvert <- OUwie.fixed(tree, trait, model=c("OUMA"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq,theta=theta, shift.point=0.5, algorithm="invert", check.identify=FALSE)
    OUMA3Point <- OUwie.fixed(tree, trait, model=c("OUMA"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq,theta=theta, shift.point=0.5, algorithm="three.point", check.identify=FALSE)
    comparison <- identical(round(as.numeric(OUMAInvert$loglik),5), round(as.numeric(OUMA3Point$loglik),5))
    expect_true(comparison)
})


test_that("testing OUMVA likelihood invert vs three.point", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    
    alpha=c(3.0793193, 0.6060786)
    sigma.sq=c(0.4735485, 1.7049102)
    theta=c(1.68189033, -1.032546)
    
    OUMVAInvert <- OUwie.fixed(tree, trait, model=c("OUMVA"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq,theta=theta, shift.point=0.5, algorithm="invert", check.identify=FALSE)
    OUMVA3Point <- OUwie.fixed(tree, trait, model=c("OUMVA"), simmap.tree=FALSE, scaleHeight=TRUE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq,theta=theta, shift.point=0.5, algorithm="three.point", check.identify=FALSE)
    comparison <- identical(round(as.numeric(OUMVAInvert$loglik),5), round(as.numeric(OUMVA3Point$loglik),5))
    expect_true(comparison)
})


test_that("testing BM1 three-point likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="BM1", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="three.point", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,5), round(-21.95911, 5))
    expect_true(comparison)
})


test_that("testing BMS three-point likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="BMS", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="three.point", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,5), round(-17.85074,5))
    expect_true(comparison)
})


test_that("testing OU1 three-point likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OU1", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="three.point", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,4), round(-21.74538,4))
    expect_true(comparison)
})


test_that("testing OUM three-point likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUM", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="three.point", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,3), round(-19.51388,3))
    expect_true(comparison)
})


test_that("testing OUMV three-point likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUMV", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="three.point", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,3), round(-14.79506,3))
    expect_true(comparison)
})


test_that("testing OUMA three-point likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUMA", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="three.point", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,3), round(-19.42678,3))
    expect_true(comparison)
})


test_that("testing OUMVA three-point likelihood", {
    skip_on_cran()
    
    data(tworegime)
    set.seed(42)
    ouwiefit <- OUwie(tree, trait, model="OUMVA", scaleHeight=TRUE, root.station=FALSE, shift.point=0.5, algorithm="three.point", quiet=TRUE, check.identify=FALSE)
    comparison <- identical(round(ouwiefit$loglik,3), round(-14.059,3))
    expect_true(comparison)
})


test_that("testing mserr vs three-point", {
    skip_on_cran()
    
    data(tworegime)
    trait[,4] <- abs(rnorm(length(tree$tip.label), mean = 1, sd = 1/10))
    alpha=c(5,9)
    sigma.sq=c(1, 2)
    theta=c(5, 10)
    INV <- c(OUwie.fixed(tree,trait,model=c("OUMVA"), simmap.tree=FALSE, scaleHeight=TRUE, tip.fog = "known", clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, algorithm = "invert", quiet=TRUE, check.identify=FALSE)$loglik)
    TPT <- OUwie.fixed(tree,trait,model=c("OUMVA"), simmap.tree=FALSE, scaleHeight=TRUE, tip.fog = "known", clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, algorithm = "three.point", quiet=TRUE, check.identify=FALSE)$loglik
    
    comparison <- identical(round(INV, 5), round(TPT, 5))
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
#    OUwie(tree,test.data,model="BM1",simmap.tree=TRUE, root.station=FALSE, algorithm="invert")
#)

    ## fit using fitContinuous
#    fitContinuous(tree,y)
#})






