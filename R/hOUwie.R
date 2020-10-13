# set of functions for the hidden rates OU model

# exported function with all the bells and whistles
hOUwie <- function(phy, data, 
                   rate.cat, rate.mat=NULL, model.cor="ARD", root.p="yang", lb=1e-3, ub=10,
                   model.ou="BM1", root.station.ou=FALSE, get.root.theta=FALSE, shift.point=0.5,
                   p=NULL, ip=NULL, nSim=1000, nCores=1){
  # check that tips and data match
  # check for invariance of tip states and not that non-invariance isn't just ambiguity
  
  #Ensures that weird root state probabilities that do not sum to 1 are input:
  if(!is.null(root.p)){
    if(!is.character(root.p)){
      root.p <- root.p/sum(root.p)
    }
  }
  
  # organize the data into the corHMM data and the OUwie data
  hOUwie.dat <- organizeHOUwieDat(data)
  nObs <- length(hOUwie.dat$ObservedTraits)
  
  #Some initial values for use later
  model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data=hOUwie.dat$data.cor,rate.cat=rate.cat, ntraits = nObs, model = model.cor)
  phy <- reorder(phy, "pruningwise")
  
  # this allows for custom rate matricies!
  if(!is.null(rate.mat)){
    order.test <- FALSE
    rate.mat[rate.mat == 0] <- NA
    rate <- rate.mat
    model.set.final$np <- max(rate, na.rm=TRUE)
    rate[is.na(rate)]=max(rate, na.rm=TRUE)+1
    model.set.final$rate <- rate
    model.set.final$index.matrix <- rate.mat
    model.set.final$Q <- matrix(0, dim(rate.mat)[1], dim(rate.mat)[2])
    ## for precursor type models ##
    col.sums <- which(colSums(rate.mat, na.rm=TRUE) == 0)
    row.sums <- which(rowSums(rate.mat, na.rm=TRUE) == 0)
    drop.states <- col.sums[which(col.sums == row.sums)]
    if(length(drop.states > 0)){
      model.set.final$liks[,drop.states] <- 0
    }
    ###############################
  }
  
  # get the appropriate OU model structure
  index.ou <- getParamStructure(model.ou, "three.point", root.station.ou, get.root.theta, model.set.final$np)
  # p is organized into 2 groups with the first set being corHMM and the second set being OUwie
  # organized as c(trans.rt, alpha, sigma.sq, theta)
  est.p <- numeric(model.set.final$np + max(index.ou, na.rm=TRUE))
  # MLE search options
  opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
  opts.quick <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.1)
  # evaluate likelihood
  if(!is.null(p)){
    cat("Calculating likelihood from a set of fixed parameters", "\n")
    out<-NULL
    est.pars<-log(p)
    out$objective <- hOUwie.dev(est.pars, phy=phy, data.cor=hOUwie.dat$data.cor , data.ou=hOUwie.dat$data.ou, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=nSim, nCores=nCores)
    loglik <- -out$objective
    est.pars <- exp(est.pars)
  }else{
    cat("This feature is not yet finzalized\n")
    out<-NULL
    start.cor <- rep(10/sum(phy$edge.length), model.set.final$np)
    alpha <- rep(1e-9, length(index.ou[1,]))
    sig2 <- rep(var(hOUwie.dat$data.ou[,3]), length(index.ou[2,]))
    theta <- rep(mean(hOUwie.dat$data.ou[,3]), length(index.ou[3,]))
    start.ou <- cbind(alpha, sig2, theta)[t(!is.na(index.ou))]
    starts = c(start.cor, start.ou)
    
    if(model.ou == "BMS" | model.ou == "BM1"){
      lower = log(c(rep(lb, model.set.final$np), rep(1e-5, length(sig2))))
      upper = log(c(rep(ub, model.set.final$np), rep(10, length(sig2))))
    }else{
      lower = log(c(rep(lb, model.set.final$np), rep(1e-9, length(alpha)), rep(1e-5, length(sig2)), rep(1e-10, length(theta))))
      upper = log(c(rep(ub, model.set.final$np), rep(10, length(alpha)), rep(10, length(sig2)), rep(10, length(theta))))
    }
    cat("Starting an initial search of parameters with a single simmap...\n")
    out = nloptr(x0=log(starts), eval_f=hOUwie.dev, lb=lower, ub=upper, opts=opts.quick, phy=phy, data.cor=hOUwie.dat$data.cor , data.ou=hOUwie.dat$data.ou, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=1, nCores=nCores)
    cat("\n\nStarting the final ML search with informed parameters...\n")
    starts <- exp(out$solution)
    out = nloptr(x0=log(starts), eval_f=hOUwie.dev, lb=lower, ub=upper, opts=opts, phy=phy, data.cor=hOUwie.dat$data.cor , data.ou=hOUwie.dat$data.ou, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=nSim, nCores=nCores)
  }
  return(out)
}

# for a single set of parameters, evaluate the hOUwie likelihood
hOUwie.dev <- function(p, phy, data.cor, data.ou, liks, Q, rate, root.p, rate.cat, index.ou, model.ou, nSim, nCores){
  
  # need to organize the parameters within this funciton. currently optimizing too many...
  # params are given in log form
  p <- exp(p)
  # define which params are for the HMM
  k <- max(rate)-1
  p.mk <- p[1:k]
  # set the OU params
  p.ou <- p[(k+1):length(p)] 
  Rate.mat <- matrix(1, 3, dim(rate)[2])
  Rate.mat[] <- c(p.ou, 1e-10)[index.ou]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  
  # fit the corHMM model. if rate.cat > 1, then ensure the order of the state mats is fast to slow.
  if(rate.cat == 1){
    Mk.loglik <- -corHMM:::dev.corhmm(p.mk, phy, liks, Q, rate, root.p, rate.cat, FALSE, FALSE)
  }else{
    Mk.loglik <- -corHMM:::dev.corhmm(p.mk, phy, liks, Q, rate, root.p, rate.cat, TRUE, FALSE)
  }
  # set up the rate matrix
  Q[] <- c(p.mk, 0)[rate]
  diag(Q) <- -rowSums(Q)
  # fit the ancestral state reconstruction
  lik.anc <- corHMM:::dev.ancRECON.marginal(p.mk, phy, liks, Q, rate, root.p, rate.cat, FALSE)
  # simulate a set of simmaps
  simmap <- mclapply(1:nSim, function(x) makeSimmap(phy, lik.anc$lik.tip.states, lik.anc$lik.anc.states, Q, 1, 1)[[1]], mc.cores = nCores)
  # fit the OU models to the simmaps
  if(model.ou == "BMS" | model.ou == "BM1"){
    OU.loglik <- mclapply(simmap, function(x) OUwie.fixed(x, data.ou, model=model.ou, simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, sigma.sq=sigma.sq, algorithm="invert", quiet = TRUE)$loglik, mc.cores = nCores)
  }else{
    OU.loglik <- mclapply(simmap, function(x) OUwie.fixed(x, data.ou, model=model.ou, simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="invert", quiet = TRUE)$loglik, mc.cores = nCores)
  }
  OU.loglik <- mean(unlist(OU.loglik))
  cat("\rpars =", round(p, 5), "lik =", OU.loglik + Mk.loglik, "                ")
  
  return(-(OU.loglik + Mk.loglik))
}

# hOUwie's input data will be similar to OUwie, with the exception that there can be more than a single trait being evaluated thus it's defined as column 1 is species name, the last column is the continuous trait and anything in between are discrete characters
organizeHOUwieDat <- function(data){
  # return a list of corHMM data and OU data
  data.cor <- data[1:dim(data)[2]-1]
  data.cor <- corHMM:::corProcessData(data.cor)
  data.ou <- data.frame(sp = data[,1], reg = data.cor$corData[,2], x = data[,3])
  return(list(StateMats = data.cor$StateMats, 
              PossibleTraits = data.cor$PossibleTraits,
              ObservedTraits = data.cor$ObservedTraits,
              data.cor = data.cor$corData,
              data.ou = data.ou))
}

# different OU models have different parameter structures. This will evaluate the appropriate one.
getParamStructure <- function(model, algorithm, root.station, get.root.theta, k){
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
      index.mat <- rbind(index.mat, rep(NA,k))
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
        max.par.so.far <- max(index.mat)
        index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
      }
      param.count <- np+k
    }
    if(root.station==FALSE){
      if(algorithm == "three.point"){
        index.mat <- rbind(index.mat, rep(NA,k))
      }
      param.count <- np+1
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
  return(index.mat)
}


# testing

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
# model.ou <- "OUMV"
# root.p <- "yang"
# 
# test <- OUwie:::hOUwie(phy, data, 1, model.ou = model.ou, ub = 3, nSim = 2)
# #debug(OUwie:::hOUwie)
# 
# p = c(0.001, 0.001, 0.1027, 0.14834, 1.33008, 0.19581, 3.52749, 4.23788)
# hOUwie.dat <- OUwie:::organizeHOUwieDat(data)
# nObs <- length(hOUwie.dat$ObservedTraits)
# model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data=hOUwie.dat$data.cor,rate.cat=rate.cat, ntraits = nObs, model = model.cor)
# phy <- reorder(phy, "pruningwise")
# index.ou <- OUwie:::getParamStructure(model.ou, "three.point", FALSE, FALSE, model.set.final$np)
# 
# OUwie:::hOUwie.dev(p = log(p), phy = phy, data.cor = OUwie:::organizeHOUwieDat(data)$data.cor, data.ou = OUwie:::organizeHOUwieDat(data)$data.ou, liks = model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat=rate.cat, index.ou=index.ou, model.ou=model.ou, nSim=1000, nCores=1)
# 
# undebug(OUwie:::hOUwie.dev)
# 
# 
# 
