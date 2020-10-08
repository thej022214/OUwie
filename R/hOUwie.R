# set of functions for the hidden rates OU model

# exported function with all the bells and whistles
hOUwie <- function(phy, data, 
                   rate.cat, rate.mat=NULL, model="ARD", root.p="yang", lb=1e-7, ub=10, 
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
  lb <- log(lb)
  ub <- log(ub)
  model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data=hOUwie.dat$data.cor,rate.cat=rate.cat, ntraits = nObs, model = model)
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
  
  # p is organized into 2 groups with the first set being corHMM and the second set being OUwie
  lower = rep(lb, model.set.final$np)
  upper = rep(ub, model.set.final$np)
  p <- c(trans.rt, alpha, sigma.sq, theta)
  
  # MLE search options
  opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
  
  # evaluate likelihood
  if(!is.null(p)){
    cat("Calculating likelihood from a set of fixed parameters", "\n")
    out<-NULL
    est.pars<-log(p)
    out$objective <- hOUwie.dev(est.pars, phy=phy, data.cor=hOUwie.dat$data.cor , data.ou=hOUwie.dat$data.ou, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat = rate.cat, nSim = nSim, nCores = nCores)
    loglik <- -out$objective
    est.pars <- exp(est.pars)
  }else{
    cat("This feature is not yet implemented\n")
    # out = nloptr(x0=log(starts), eval_f=dev.corhmm, lb=lower, ub=upper, opts=opts, phy=phy, liks=model.set.final$liks,Q=model.set.final$Q,rate=model.set.final$rate,root.p=root.p, rate.cat = rate.cat, order.test = order.test, lewis.asc.bias = lewis.asc.bias)
  }
  return(list(loglik = loglik, est.pars = est.pars))
}

# for a single set of parameters, evaluate the hOUwie likelihood
hOUwie.dev <- function(p, phy, data.cor, data.ou, liks, Q, rate, root.p, rate.cat, nSim, nCores){
  
  # params are given in log form
  p <- exp(p)
  # define which params are for the HMM
  p.mk = p[1:(max(rate)-1)]
  # set the OU params
  alpha = p[3:4]
  sigma.sq = p[5:6]
  theta = p[7:8]
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
  OU.loglik <- mclapply(simmap, function(x) OUwie.fixed(x, data.ou, model=c("OUMVA"), simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="invert")$loglik, mc.cores = nCores)
  OU.loglik <- mean(unlist(OU.loglik))
  
  return(OU.loglik + Mk.loglik)
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
# model <- "ARD"
# root.p <- "yang"
# 
# trans.rt=c(0.01607089, 0.01707089)
# alpha=c(0.5632459,0.1726052)
# sigma.sq=c(0.1064417,0.3461386)
# theta=c(1.678196,0.4185894)
# 
# p <- c(trans.rt, alpha, sigma.sq, theta)
# hOUwie.dat <- organizeHOUwieDat(data)
# nObs <- length(hOUwie.dat$ObservedTraits)
# lb <- log(1e-7)
# ub <- log(10)
# model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy=phy,data=hOUwie.dat$data.cor,rate.cat=rate.cat, ntraits = nObs, model = model)
# phy <- reorder(phy, "pruningwise")
# est.pars<-log(p)
# 
# nSim <- 100
# nCores <- 2
# 
# hOUwie.dev(est.pars, phy=phy, data.cor=hOUwie.dat$data.cor , data.ou=hOUwie.dat$data.ou, liks=model.set.final$liks, Q=model.set.final$Q, rate=model.set.final$rate, root.p=root.p, rate.cat = rate.cat, nSim = nSim, nCores = nCores)

