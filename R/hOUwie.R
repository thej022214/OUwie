# set of functions for the hidden rates OU model

# exported function with all the bells and whistles
#'@author James Boyko
#'@param phy a phylogeny of class phylo
#'@param data data.frame where [,1] = species name,  [,2:j] are discrete traits, [,j:m] are the continuous trait (and possible mserr). if mserr = "none" j=m-1, if mserr = "known j=m-2. m is the number of total columns in the dataset.
#'@param rate.cat the number rate categories in the corHMM model
#'@param model.cor a single corHMM model or a vector of corHMM models (ARD, ER, or SYM) 
#'@param index.cor custom rate/index matrix
#'@param lb.cor lower bound for trans rate
#'@param ub.cor upper bound for trans rate
#'@param root.p designates the root prior (yang, maddfitz, flat, numeric vector)
#'@param model.ou designates the default OU model (BM1, BMS, OU1, OUM, OUMV, OUMA, OUMVA)
#'@param index.ou designates which ou params are being estimated
#'@param root.station logical indicating whether to assume a random starting point (TRUE) or a fixed starting point (FALSE) 
#'@param get.root.theta logical indicating whether the starting state, theta_0, should be estimated
#'@param mserr designates whether mserr is included (known) or not (none)
#'@param lb.ou designates the lower bounds for alpha, sigma.sq, optima
#'@param ub.ou designates the upper bounds for alpha, sigma.sq, optima
#'@param p a vector of params to calculate a fixed likelihood
#'@param ip a vector of initial params 
#'@param nSim the number of simmaps a single set of params is evaluated over
#'@param opts optioins for nloptr
#'@param nCores number of parallel cores
#'@param weighted a logicial indicating whether to take the average OU likelihood over all simmaps (FALSE) or just take the maximum likelihood
#'@param quiet a logical indicating whether to output user messages
hOUwie <- function(phy, data, rate.cat, 
                   model.cor="ER", index.cor=NULL, root.p="yang", lb.cor=1e-5, ub.cor=10,
                   model.ou="BM1", index.ou=NULL, root.station=FALSE, get.root.theta=FALSE, mserr = "none", lb.ou=c(1e-5,1e-5,1e-5), ub.ou=c(21,21,21),
                   p=NULL, ip=NULL, nSim=1000, opts=NULL, nCores=1, weighted=FALSE, quiet=FALSE){
  # check that tips and data match
  # check for invariance of tip states and not that non-invariance isn't just ambiguity
  
  #Ensures that weird root state probabilities that do not sum to 1 are input:
  if(!is.null(root.p)){
    if(!is.character(root.p)){
      root.p <- root.p/sum(root.p)
    }
  }
  algorithm="three.point"
  # organize the data into the corHMM data and the OUwie data
  # need to add a way to shift negative continuous variables to positive then shift back
  hOUwie.dat <- organizeHOUwieDat(data, mserr)
  nObs <- length(hOUwie.dat$ObservedTraits)
  # a way to speed up the three.point function
  tip.paths <- lapply(1:length(phy$tip.label), function(x) OUwie:::getPathToRoot(phy, x))
  Tmax <- max(branching.times(phy))
  
  #scale the tree to a root height of 1
  # phy$edge.length <- phy$edge.length/max(branching.times(phy))
  
  model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data=hOUwie.dat$data.cor,rate.cat=rate.cat, ntraits = nObs, model = model.cor)
  # get the appropriate OU model structure
  if(is.null(index.ou)){
    index.ou <- getOUParamStructure(model.ou, "three.point", root.station, get.root.theta, dim(model.set.final$Q)[1])
  }
  #reorder phy
  phy <- reorder(phy, "pruningwise")
  
  # this allows for custom rate matricies!
  if(!is.null(index.cor)){
    order.test <- FALSE
    index.cor[index.cor == 0] <- NA
    rate <- index.cor
    model.set.final$np <- max(rate, na.rm=TRUE)
    rate[is.na(rate)]=max(rate, na.rm=TRUE)+1
    model.set.final$rate <- rate
    model.set.final$index.matrix <- index.cor
    model.set.final$Q <- matrix(0, dim(index.cor)[1], dim(index.cor)[2])
    ## for precursor type models ##
    col.sums <- which(colSums(index.cor, na.rm=TRUE) == 0)
    row.sums <- which(rowSums(index.cor, na.rm=TRUE) == 0)
    drop.states <- col.sums[which(col.sums == row.sums)]
    if(length(drop.states > 0)){
      model.set.final$liks[,drop.states] <- 0
    }
    ## need to do anything to the ouwie matrix?
    ###############################
  }

  # default MLE search options
  if(is.null(opts)){
    opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
  }
  # p is organized into 2 groups with the first set being corHMM and the second set being OUwie
  # organized as c(trans.rt, alpha, sigma.sq, theta)
  # evaluate likelihood
  if(!is.null(p)){
    cat("Calculating likelihood from a set of fixed parameters", "\n")
    out<-NULL
    est.pars<-log(p)
    out$objective <- hOUwie.dev(est.pars, 
                                phy=phy, rate.cat=rate.cat,
                                data.cor=hOUwie.dat$data.cor, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, 
                                data.ou=hOUwie.dat$data.ou, index.ou=index.ou, algorithm=algorithm, mserr=mserr,
                                nSim=nSim, nCores=nCores, tip.paths=tip.paths, weighted=weighted)
    loglik <- -out$objective
    est.pars <- exp(est.pars)
  }else{
    cat("This feature is not yet finzalized\n")
    out<-NULL
    means.by.regime <- with(hOUwie.dat$data.ou, tapply(hOUwie.dat$data.ou[,3], hOUwie.dat$data.ou[,2], mean))
    if(length(unique(na.omit(index.ou[3,]))) > length(means.by.regime)){
      start.theta <- rep(means.by.regime, length(unique(index.ou[3,]))/length(means.by.regime))
    }else{
      start.theta <- rep(mean(hOUwie.dat$data.ou[,3]), length(unique(na.omit(index.ou[3,]))))
    }
    start.cor <- rep(10/sum(phy$edge.length), model.set.final$np)
    start.ou <- c(rep(log(2)/Tmax, length(unique(na.omit(index.ou[1,])))), 
                  rep(var(hOUwie.dat$data.ou[,3]), length(unique(na.omit(index.ou[2,])))), 
                  start.theta)
    starts = c(start.cor, start.ou)
    lower = log(c(rep(lb.cor, model.set.final$np), 
                  rep(lb.ou[1], length(unique(na.omit(index.ou[1,])))), 
                  rep(lb.ou[2], length(unique(na.omit(index.ou[2,])))), 
                  rep(lb.ou[3], length(unique(na.omit(index.ou[3,]))))))
    upper = log(c(rep(ub.cor, model.set.final$np), 
                  rep(ub.ou[1], length(unique(na.omit(index.ou[1,])))), 
                  rep(ub.ou[2], length(unique(na.omit(index.ou[2,])))), 
                  rep(ub.ou[3], length(unique(na.omit(index.ou[3,]))))))
    cat("\nStarting a search of parameters with", nSim, "simmaps...\n")
    out = nloptr(x0=log(starts), eval_f=hOUwie.dev, lb=lower, ub=upper, opts=opts, 
                 phy=phy, rate.cat=rate.cat,
                 data.cor=hOUwie.dat$data.cor, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, 
                 data.ou=hOUwie.dat$data.ou, index.ou=index.ou, algorithm=algorithm, mserr=mserr,
                 nSim=nSim, nCores=nCores, tip.paths=tip.paths, weighted=weighted)
    cat("\nFinished.\n")
  }
  return(out)
}

# for a single set of parameters, evaluate the hOUwie likelihood
hOUwie.dev <- function(p, phy, rate.cat, 
                       data.cor, liks, Q, rate, root.p,
                       data.ou, index.ou, algorithm, mserr,
                       nSim, nCores, tip.paths=NULL, weighted = FALSE){
  # params are given in log form
  p <- exp(p)
  # define which params are for the HMM
  k <- max(rate)-1
  p.mk <- p[1:k]
  # set the OU params
  p.ou <- p[(k+1):length(p)] 
  Rate.mat <- matrix(1, 3, dim(rate)[2])
  index.ou[is.na(index.ou)] <- max(index.ou, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[index.ou]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  if(algorithm == "invert"){
    theta = NULL
  }
  # fit the corHMM model. if rate.cat > 1, then ensure the order of the state mats is fast to slow.
  Mk.loglik <- -corHMM:::dev.corhmm(log(p.mk), phy, liks, Q, rate, root.p, rate.cat, TRUE, FALSE)
  # set up the rate matrix
  Q[] <- c(p.mk, 0)[rate]
  diag(Q) <- -rowSums(Q)
  # fit the ancestral state reconstruction
  # simulate a set of simmaps
  simmap <- corHMM:::makeSimmap(phy, data.cor, Q, rate.cat, nSim = nSim)
  # fit the OU models to the simmaps
  
  OU.loglik <- mclapply(simmap, function(x) OUwie.basic(x, data.ou, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm=algorithm, tip.paths=tip.paths, mserr=mserr), mc.cores = nCores)
  
  if(weighted == TRUE){
    # get the likelihoods of the simmaps
    OU.loglik <- max(unlist(OU.loglik))
    # OU.loglik <- log(mean(exp(unlist(OU.loglik)-comp)))+comp
    cat("\rpars =", round(p, 5), "lik =", -(OU.loglik + Mk.loglik), "                ")
    return(-(OU.loglik + Mk.loglik))
  }else{
    # comp <- max(unlist(OU.loglik))
    OU.loglik <- log(mean(exp(unlist(OU.loglik))))
    # OU.loglik <- log(mean(exp(unlist(OU.loglik)-comp)))+comp
    cat("\rpars =", round(p, 5), "& lik =", -(OU.loglik + Mk.loglik), "                ")
    return(-(OU.loglik + Mk.loglik))
  }
}

# hOUwie's input data will be similar to OUwie, with the exception that there can be more than a single trait being evaluated thus it's defined as column 1 is species name, the last column is the continuous trait and anything in between are discrete characters
organizeHOUwieDat <- function(data, mserr){
  # return a list of corHMM data and OU data
  if(mserr=="known"){
    data.cor <- data[, 1:(dim(data)[2]-2)]
    data.cor <- corHMM:::corProcessData(data.cor)
    data.ou <- data.frame(sp = data[,1], 
                          reg = data.cor$corData[,2], 
                          x = data[, dim(data)[2]-1],
                          err = data[, dim(data)[2]])
  }
  if(mserr=="none"){
    data.cor <- data[, 1:(dim(data)[2]-1)]
    data.cor <- corHMM:::corProcessData(data.cor)
    data.ou <- data.frame(sp = data[,1], 
                          reg = data.cor$corData[,2], 
                          x = data[, dim(data)[2]])
  }
  return(list(StateMats = data.cor$StateMats, 
              PossibleTraits = data.cor$PossibleTraits,
              ObservedTraits = data.cor$ObservedTraits,
              data.cor = data.cor$corData,
              data.ou = data.ou))
}

# different OU models have different parameter structures. This will evaluate the appropriate one.
getOUParamStructure <- function(model, algorithm, root.station, get.root.theta, k){
  index.mat <- matrix(0, 2,k)
  if(algorithm == "three.point"){
    Rate.mat <- matrix(1, 3, k)
  }else{
    Rate.mat <- matrix(1, 2, k)
  }
  if (model == "BM1"){
    np <- 1
    index.mat[1,1:k] <- NA
    index.mat[2,1:k] <- 1
    if(algorithm == "three.point"){
      index.mat <- rbind(index.mat, rep(2,k))
    }
    param.count <- np+1
  }
  #The group mean model of Thomas et al, is trivial: set root.station to be TRUE:
  if (model == "BMS"){
    np <- k
    index.mat[1,1:k] <- NA
    index.mat[2,1:k] <- 1:np
    if(root.station==TRUE){
      if(algorithm == "three.point"){
        max.par.so.far <- max(index.mat, na.rm = TRUE)
        index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
      }
      param.count <- np+k
    }
    if(root.station==FALSE){
      if(algorithm == "three.point"){
        max.par.so.far <- max(index.mat, na.rm = TRUE)
        index.mat <- rbind(index.mat, max.par.so.far+1)
      }
      param.count <- np+k
    }
  }
  if (model == "OU1"){
    np <- 2
    index.mat[1,1:k] <- 1
    index.mat[2,1:k] <- 2
    if(algorithm == "three.point"){
      max.par.so.far <- max(index.mat)
      index.mat <- rbind(index.mat, rep(3,k))
    }
    param.count <- np + 1
    if(get.root.theta == TRUE){
      param.count <- param.count + 1
    }
  }
  if (model == "OUM"){
    np <- 2
    index.mat[1,1:k] <- 1
    index.mat[2,1:k] <- 2
    if(algorithm == "three.point"){
      max.par.so.far <- max(index.mat)
      index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
    }
    param.count <- np + k
    if(get.root.theta == TRUE){
      param.count <- param.count + 1
    }
  }
  if (model == "OUMV") {
    np <- k+1
    index.mat[1,1:k] <- 1
    index.mat[2,1:k] <- 2:(k+1)
    if(algorithm == "three.point"){
      max.par.so.far <- max(index.mat)
      index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
    }
    param.count <- np + k
    if(get.root.theta == TRUE){
      param.count <- param.count + 1
    }
  }
  if (model == "OUMA") {
    np <- k+1
    index.mat[1,1:k] <- 1:k
    index.mat[2,1:k] <- k+1
    if(algorithm == "three.point"){
      max.par.so.far <- max(index.mat)
      index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
    }
    param.count <- np+k
    if(get.root.theta == TRUE){
      param.count <- param.count + 1
    }
  }
  if (model == "OUMVA") {
    np <- k*2
    index <- matrix(TRUE,2,k)
    index.mat[index] <- 1:(k*2)
    if(algorithm == "three.point"){
      max.par.so.far <- max(index.mat)
      index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
    }
    param.count <- np+k
    if(get.root.theta == TRUE){
      param.count <- param.count + 1
    }
  }
  if (model == "TrendyM"){
    index.mat <- matrix(0,3,k)
    np <- k+2 #We have k regimes, but 1 sigma, plus the root value
    index.mat[1,1:k] <- np+1
    index.mat[2,1:k] <- 1
    index.mat[3,1:k] <- 2:(np-1)
    param.count<-np
    root.station=TRUE
    trendy=1
  }
  if (model == "TrendyMS"){
    index.mat <- matrix(0,3,k)
    np <- (k*2)+1 #We have k regimes and assume k trends and k sigmas, plus the root value
    index.mat[1,1:k] <- np+1
    index.mat[2,1:k] <- 1:k
    index.mat[3,1:k] <- (k+1):(np-1)
    param.count <- np
    root.station=FALSE
    trendy=1
  }
  if(algorithm=="three.point"){
    rownames(index.mat) <- c("alpha", "sigma2", "theta")
  }
  return(index.mat)
}

# a basic version of OUwie that will evaluate params based on an index matrix
OUwie.basic<-function(phy, data, simmap.tree=TRUE, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, get.root.theta=FALSE, shift.point=0.5, alpha, sigma.sq, theta, mserr="none", algorithm="three.point", tip.paths=NULL){
  
  # organize tip states based on what the simmap suggests
  mapping <- unlist(lapply(phy$maps, function(x) names(x[length(x)])))
  nTip <- length(phy$tip.label)
  TipStates <- mapping[match(match(data[,1], phy$tip.label), phy$edge[,2])]
  data[,2] <- TipStates
  
  #Makes sure the data is in the same order as the tip labels
  if(mserr=="none"){
    data <- data.frame(data[,2], data[,3], row.names=data[,1])
    data <- data[phy$tip.label,]
  }
  if(mserr=="known"){
    # algorithm = "invert"
    if(!dim(data)[2]==4){
      stop("You specified measurement error should be incorporated, but this information is missing.", call. = FALSE)
    }
    else{
      if(is.factor(data[,4]) == TRUE){
        stop("Check the format of the measurement error column. It's reading as a factor.",  call. = FALSE)
      }else{
        data <- data.frame(data[,2], data[,3], data[,4], row.names=data[,1])
        data <- data[phy$tip.label,]
      }
    }
  }
  
  #Values to be used throughout
  n <- max(phy$edge[,1])
  ntips <- length(phy$tip.label)

  # setup values when simmap (always simmap for hOUwie)
  k <- length(colnames(phy$mapped.edge))
  tot.states <- factor(colnames(phy$mapped.edge))
  tip.states <- factor(data[,1])
  data[,1] <- as.numeric(tip.states)
  #Obtains the state at the root
  root.state <- which(colnames(phy$mapped.edge)==names(phy$maps[[1]][1]))
  ##Begins the construction of the edges matrix -- similar to the ouch format##
  edges <- cbind(c(1:(n-1)),phy$edge,OUwie:::MakeAgeTable(phy, root.age=root.age))
  if(scaleHeight == TRUE){
    Tmax <- max(OUwie:::MakeAgeTable(phy, root.age=root.age))
    edges[,4:5]<-edges[,4:5]/Tmax
    root.age <-  1
    phy$maps <- lapply(phy$maps, function(x) x/Tmax)
  }
  edges <- edges[sort.list(edges[,3]),]
  #Resort the edge matrix so that it looks like the original matrix order
  edges <- edges[sort.list(edges[,1]),]

  if(algorithm == "three.point"){
    x <- data[,2]
    names(x) <- rownames(data)
  }else{
    x <- as.matrix(data[,2])
  }
  
  if(scaleHeight==TRUE){
    phy$edge.length <- phy$edge.length/Tmax
    Tmax <- 1
  }
  
  if(algorithm == "three.point"){
    if(simmap.tree == FALSE){
      map <- getMapFromNode(phy, tip.states, node.states, shift.point)
    }else{
      map <- phy$maps
    }
  }
  
  Rate.mat <- rbind(alpha, sigma.sq, theta)
  
  if(get.root.theta == TRUE){
    W <- OUwie:::weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=FALSE, shift.point=shift.point)
  }else{
    W <- OUwie:::weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=TRUE, shift.point=shift.point)
  }
  
  #Likelihood function for estimating model parameters
  pars <- matrix(c(theta, sigma.sq, alpha), length(theta), 3, dimnames = list(1:length(sigma.sq), c("opt", "sig", "alp")))
  if(get.root.theta == TRUE){
    expected.vals <- colSums(t(W) * c(theta0, pars[,1]))
    names(expected.vals) <- phy$tip.label
  }else{
    expected.vals <- colSums(t(W) * pars[,1])
    names(expected.vals) <- phy$tip.label
  }
  transformed.tree <- OUwie:::transformPhy(phy, map, pars, tip.paths)
  # generate a map from node based reconstructions
  if(mserr=="known"){
    TIPS <- transformed.tree$tree$edge[,2] <= length(transformed.tree$tree$tip.label)
    transformed.tree$tree$edge.length[TIPS] <- transformed.tree$tree$edge.length[TIPS] + (data[,3]^2/transformed.tree$diag/transformed.tree$diag)
  }
  try(comp <- phylolm:::three.point.compute(transformed.tree$tree, x, expected.vals, transformed.tree$diag), silent=TRUE)
  if(is.na(comp[1])){
    return(10000000)
  }else{
    neg.logl <- -as.numeric(Ntip(phy) * log(2 * pi) + comp$logd + (comp$PP - 2 * comp$QP + comp$QQ))/2
    return(neg.logl)
  }
}


# simulate a hOUwie model
hOUwie.sim <- function(phy, Q, root.freqs, alpha, sig2, theta0, theta, nMap=1){
  # simulate an Mk dataset
  dat.cor <- rTraitDisc(phy, Q, states = 1:dim(Q)[1], root.value = sample(1:dim(Q)[1], 1, prob = root.freqs))
  while(!all(levels(dat.cor) %in% dat.cor)){
    dat.cor <- rTraitDisc(phy, Q, states = 1:dim(Q)[1], root.value = sample(1:dim(Q)[1], 1, prob = root.freqs))
  }
  dat.cor <- data.frame(sp=names(dat.cor), d=dat.cor)
  # simulate a stochastic map with true Q
  simmap <- corHMM:::makeSimmap(phy, dat.cor, Q, 1, nSim=nMap)
  # lik <- corHMM:::getSimmapLik(simmap, Q)
  # simulate the ou dataset
  dat.ou <- lapply(simmap, function(x) OUwie.sim(x, simmap.tree = TRUE, alpha = alpha, sigma.sq = sig2, theta0 = theta0, theta = theta)[,2])
  dat.ou <- colMeans(do.call(rbind, dat.ou))
  # dat.ou <- OUwie.sim(simmap, simmap.tree = TRUE, alpha = alpha, sigma.sq = sig2, theta0 = theta0, theta = theta)
  # return true params and data
  data <- data.frame(sp = phy$tip.label, reg = dat.cor[,2], x = dat.ou)
  return(list(data = data, simmap = simmap))
}

# testing
# 
# require(corHMM)
# require(OUwie)
# require(parallel)
# data(tworegime)
# 
# tree$node.label <- NULL
# data <- trait
# phy <- tree
# rate.cat <- 1
# model.cor <- "ARD"
# model.ou <- "OUM"
# root.p <- "yang"
# 
# test <- OUwie:::hOUwie(phy, data, 1, model.ou = model.ou, ub = 3, nSim = 10, weighted = FALSE)
# undebug(OUwie:::hOUwie.dev)
# 
# OUwie:::hOUwie.dev(p = log(p), phy = phy, data.cor = OUwie:::organizeHOUwieDat(data)$data.cor, data.ou = OUwie:::organizeHOUwieDat(data)$data.ou, liks = model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=100, nCores=1)

# hOUwie.dev(p = log(p), phy = phy, data.cor = OUwie:::organizeHOUwieDat(data)$data.cor, data.ou = OUwie:::organizeHOUwieDat(data)$data.ou, liks = model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=100, nCores=1)


# undebug(OUwie:::hOUwie.dev)
#
# Rprof(filename = "~/2020_hOUwie/Rprof.out", append = FALSE, line.profiling = TRUE)
# summaryRprof("2020_hOUwie/Rprof.out")
#
# simulation tests
# require(corHMM)
# require(OUwie)
# require(parallel)
# data(tworegime)
#
# data <- trait
# phy <- tree
# phy$node.label <- NULL
# root.p = c(0.5, 0.5)
# p.mk <- c(0.5, 0.5)
# alpha = c(5, 5)
# sig2= c(0.1, 0.1)
# theta = c(3, 8)
# Q = matrix(c(-p.mk[1],p.mk[2],p.mk[1],-p.mk[2]), 2, 2)
# theta0 = 5
#
# data <- OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta)[[2]]
# test <- OUwie:::hOUwie(phy, data, 1, model.ou = "OUM", ub = 3, nSim = 10, weighted = TRUE)
# test <- OUwie:::hOUwie(phy, data, 1, model.ou = "OUM", ub = 3, nSim = 10, weighted = FALSE)
#
#
# OUwie:::hOUwie(phy, data, 1, model.ou = "OUM", p = c(p.mk, alpha[1], sig2[1], theta), nSim = 1000)

#
# undebug(OUwie:::hOUwie.dev)

#
# require(corHMM)
# require(OUwie)
# require(parallel)
# require(geiger)
# require(proftools)
# 
# phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = 250)
# phy <- drop.extinct(phy)
# root.p = c(0.5, 0.5)
# p.mk <- c(0.1, 0.1)
# alpha = c(1, 2)
# sig2= c(0.1, 0.1)
# theta = c(3, 8)
# Q = matrix(c(-p.mk[1],p.mk[2],p.mk[1],-p.mk[2]), 2, 2)
# theta0 = 5
# rate.cat = 1
# model.cor = "ER"
# model.ou = "OUMA"
# 
# data <- OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta)[[1]]
# p = c(0.1, 0.01, 0.1, 1, 3, 8)
# hOUwie.dat <- OUwie:::organizeHOUwieDat(data)
# nObs <- length(hOUwie.dat$ObservedTraits)
# model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data=hOUwie.dat$data.cor,rate.cat=rate.cat, ntraits = nObs, model = model.cor)
# phy <- reorder(phy, "pruningwise")
# index.ou <- OUwie:::getParamStructure(model.ou, "three.point", FALSE, FALSE, dim(model.set.final$Q)[2])
# # phy$edge.length <- phy$edge.length/max(branching.times(phy))
# 
# OUwie:::hOUwie.dev(p = log(p), phy = phy, data.cor = OUwie:::organizeHOUwieDat(data)$data.cor, data.ou = OUwie:::organizeHOUwieDat(data)$data.ou, liks = model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=10, nCores=1, algorithm = "three.point")

# # # weigthed <- sapply(1:100, function(x) OUwie:::hOUwie.dev(p = log(p), phy = phy, data.cor = OUwie:::organizeHOUwieDat(data)$data.cor, data.ou = OUwie:::organizeHOUwieDat(data)$data.ou, liks = model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=10, nCores=1, weighted = TRUE))
# # 
# # # unweigthed <- sapply(1:100, function(x) OUwie:::hOUwie.dev(p = log(p), phy = phy, data.cor = OUwie:::organizeHOUwieDat(data)$data.cor, data.ou = OUwie:::organizeHOUwieDat(data)$data.ou, liks = model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=100, nCores=1))
# # 
# # 
# # pd1 <- profileExpr(OUwie:::hOUwie.dev(p = log(p), phy = phy, data.cor = OUwie:::organizeHOUwieDat(data)$data.cor, data.ou = OUwie:::organizeHOUwieDat(data)$data.ou, liks = model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=100, nCores=1, algorithm = "three.point"))
# # 
# # pd4 <- profileExpr(OUwie:::hOUwie.dev(p = log(p), phy = phy, data.cor = OUwie:::organizeHOUwieDat(data)$data.cor, data.ou = OUwie:::organizeHOUwieDat(data)$data.ou, liks = model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=100, nCores=4, algorithm = "three.point"))
# # 
# # pd10 <- profileExpr(OUwie:::hOUwie.dev(p = log(p), phy = phy, data.cor = OUwie:::organizeHOUwieDat(data)$data.cor, data.ou = OUwie:::organizeHOUwieDat(data)$data.ou, liks = model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=100, nCores=10, algorithm = "three.point"))
# # 
# # 
# # 
# # hotPaths(pd1, total.pct = 10.0)
# # hotPaths(pd4, total.pct = 10.0)
# # hotPaths(pd10, total.pct = 10.0)
# # 
# # # weighted.search <- OUwie:::hOUwie(phy, data, rate.cat, model.cor = model.cor, model.ou = model.ou, weighted = TRUE, nSim = 10)
# # unweighted.search <- OUwie:::hOUwie(phy, data, rate.cat, model.cor = model.cor, model.ou = model.ou, weighted = FALSE, nSim = 100)