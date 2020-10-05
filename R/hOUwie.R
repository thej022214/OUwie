# set of functions for the hidden rates OU model

# exported function with all the bells and whistles
hOUwie <- function(phy, data){
  
  # organize the data into the corHMM data and the OUwie data
  hOUwie.dat <- organizeHOUwieDat(trait)
  
  # p is organized into 2 groups with the first set being corHMM and the second set being OUwie
  p <- c(trans.rt, alpha, sigma.sq, theta)
  
  # evaluate likelihood
  # repeat
  
}

# for a single set of parameters, evaluate the hOUwie likelihood
hOUwie.dev <- function(p, phy, data.cor, data.ou, model, nSim){
  
  # define which params are for the HMM
  p.mk = p[1:2]
  # fit the corHMM model
  Mk <- corHMM(phy = tree, data = data.cor, rate.cat = 1, model = model, p = p.mk, get.tip.states = TRUE)
  MK$phy$node.label <- NULL
  
  # simulate a set of simmaps
  simmap <- makeSimmap(MK$phy, MK$tip.states, MK$states, trans.mat, nSim, nCores=1)
  
  # set the OU params
  alpha = p[3:4]
  sigma.sq = p[5:6]
  theta = p[7:8]
  # fit the OU models to the simmaps
  OURes <- lapply(simmap, function(x) OUwie.fixed(x,trait,model=c("OUMVA"), simmap.tree=TRUE, scaleHeight=FALSE, clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, algorithm="invert"))
  
  # combine the likelihoods
  OU.loglik <- unlist(lapply(OURes, function(x) x$loglik))
  Mk.loglik <- MK$loglik
  
  return(mean(OU.loglik) + Mk.loglik)
}

# hOUwie's input data will be similar to OUwie, with the exception that there can be more than a single trait being evaluated thus it's defined as column 1 is species name, the last column is the continuous trait and anything in between are discrete characters
organizeHOUwieDat <- function(data){
  # return a list of corHMM data and OU data
  data.cor <- data[1:dim(data)[2]-1]
  data.ou <- data
  return(list(data.cor, data.ou))
}


# testing
# require(corHMM)
# require(OUwie)
# require(parallel)
# data(tworegime)
# 
# trans.rt=c(0.01607089, 0.01707089)
# alpha=c(0.5632459,0.1726052)
# sigma.sq=c(0.1064417,0.3461386)
# theta=c(1.678196,0.4185894)
# 
# hOUwie.dat <- organizeHOUwieDat(data = trait)
# p <- c(trans.rt, alpha, sigma.sq, theta)
# 
# hOUwie.dev(p, tree, hOUwie.dat[[1]], hOUwie.dat[[2]], "ARD", 20)
# 



