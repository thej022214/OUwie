context("test-houwie-sampling.R")

## The map sampling machinery is where hOUwie decides which branches and which
## observations belong together. A mistake there does not raise an error, it just
## reweights the wrong branch or pairs a likelihood with a map that did not produce it.
## These check two invariants that would otherwise go unnoticed: the answer cannot
## depend on the order the data were typed in, and two representations of the same
## root history must give the same expectations.

make.sampling.tree <- function(){
    set.seed(42)
    phy <- rcoal(8)
    phy$edge.length <- phy$edge.length / max(branching.times(phy))
    phy
}

make.sampling.data <- function(phy){
    data.frame(sp = phy$tip.label,
               reg = c(1, 1, 2, 2, 1, 2, 1, 2),
               x = c(1.2, 0.9, 3.1, 3.4, 1.0, 2.9, 1.1, 3.2))
}

## a two segment root edge and a one segment root edge that spend the same time in the
## same regime describe the same history, so they have to give the same expectations
make.root.edge.simmap <- function(phy, split.root){
    maps <- vector("list", nrow(phy$edge))
    maps[[1]] <- if(split.root) setNames(c(0.4, 0.6), c("1", "1")) else setNames(1, "1")
    maps[[2]] <- setNames(1, "1")
    maps[[3]] <- setNames(1, "2")
    maps[[4]] <- setNames(2, "1")
    simmap <- phy
    simmap$maps <- maps
    mapped.edge <- t(sapply(maps, function(x)
        c(sum(x[names(x) == "1"]), sum(x[names(x) == "2"]))))
    colnames(mapped.edge) <- c("1", "2")
    simmap$mapped.edge <- mapped.edge
    class(simmap) <- c("simmap", "phylo")
    simmap
}

test_that("tip sampling weights each branch with its own species", {
    skip_on_cran()

    phy <- make.sampling.tree()
    houwie.dat <- make.sampling.data(phy)
    p <- c(0.5, 0.5, 0.3, 0.3, 1.0, 1.0, 1.5, 3.0, 2.0, 2.5)
    evaluate <- function(dat){
        set.seed(7)
        hOUwie(phy, dat, rate.cat = 2, discrete_model = "ER",
               continuous_model = "OUM", nSim = 10, p = p,
               sample_tips = TRUE, quiet = TRUE)$loglik
    }

    set.seed(1)
    shuffled <- houwie.dat[sample(nrow(houwie.dat)),]

    expect_equal(evaluate(shuffled), evaluate(houwie.dat))
})

test_that("a root optimum is taken from the rootward end of the root edge", {
    skip_on_cran()

    phy <- read.tree(text = "((t1:1,t2:1):1,t3:2);")
    dat <- data.frame(sp = c("t1", "t2", "t3"), reg = c(1, 2, 1),
                      x = c(1.1, 2.7, 1.4))
    expectations <- function(split.root){
        OUwie.basic(make.root.edge.simmap(phy, split.root), dat,
                    simmap.tree = TRUE, scaleHeight = FALSE, get.root.theta = TRUE,
                    alpha = c(1.2, 1.2), sigma.sq = c(0.5, 0.5), theta = c(1.0, 3.0),
                    algorithm = "three.point", return.expected.vals = TRUE)
    }

    split.root <- expectations(TRUE)

    expect_equal(length(split.root), Ntip(phy))
    expect_equal(split.root, expectations(FALSE))
})
