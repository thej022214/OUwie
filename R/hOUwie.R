# set of functions for the hidden rates OU model
##### Main exported functions ##### 
hOUwie <- function(phy, data, rate.cat, discrete_model, continuous_model, null.model=FALSE, nSim=100, root.p="yang", dual = FALSE, collapse = TRUE, root.station=FALSE, get.root.theta=FALSE, mserr = "none", lb_discrete_model=NULL, ub_discrete_model=NULL, lb_continuous_model=NULL, ub_continuous_model=NULL, recon=FALSE, nodes="internal", p=NULL, ip=NULL, optimizer="nlopt_ln", opts=NULL, quiet=FALSE, sample_tips=FALSE, sample_nodes=TRUE, adaptive_sampling=TRUE, diagn_msg=FALSE, n_starts = 1, ncores = 1){
  start_time <- Sys.time()
  # if the data has negative values, shift it right - we will shift it back later
  negative_values <- FALSE
  if(any(phy$edge.length < 0)){
    stop("Your phylogeny has negative edge lengths. I don't know what can cause this, but I know it's not good.")
  }
  if(!is.binary(phy)){
    phy <- multi2di(phy)
    warning("Your phylogeny is not bifurcating, forcing a binary tree with ape:::multi2di.")
  }
  if(any(phy$edge.length == 0)){
    phy$edge.length[phy$edge.length == 0] <- 1e-5
    warning("Your phylogeny edge lengths of 0. Adding 1e-5")
  }
  if(mserr == "none"){
    cor_dat <- data[,c(1:(dim(data)[2]-1))]
    if(any(data[,dim(data)[2]] < 0)){
      if(!quiet){
        cat("Negative values detected... adding 50 to the trait mean for optimization purposes\n")
      }
      negative_values <- TRUE
      data[,dim(data)[2]] <- data[,dim(data)[2]] + 50
    }
  }else{
    cor_dat <- data[,c(1:(dim(data)[2]-2))]
    if(any(data[,dim(data)[2]-1] < 0)){
      if(!quiet){
        cat("Negative values detected... adding 50 to the trait mean for optimization purposes\n")
      }
      negative_values <- TRUE
      data[,dim(data)[2]-1] <- data[,dim(data)[2]-1] + 50
    }
  }
  # check that tips and data match
  # check for invariance of tip states and not that non-invariance isn't just ambiguity
  if(!is.null(phy$node.label)){
    if(!quiet){
      cat("Your phylogeny had node labels, these have been removed.\n")
    }
    phy$node.label <- NULL
  }
  
  
  if(ncores > n_starts){
    cat("You have specified more cores are to be used than the number of starts. Setting ncores to be equal to the number of optimizations.\n")
    ncores <- n_starts
  }
  
  # organize the data
  phy <- reorder.phylo(phy, "pruningwise")
  hOUwie.dat <- organizeHOUwieDat(data, mserr, collapse)
  nStates <- as.numeric(max(hOUwie.dat$data.cor[,2]))
  nCol <- dim(data)[2] - ifelse(mserr == "none", 2, 3)
  Tmax <- max(branching.times(phy))
  all.paths <- lapply(1:(Nnode(phy) + Ntip(phy)), function(x) getPathToRoot(phy, x))
  
  if(class(discrete_model)[1] == "character"){
    index.disc <- getDiscreteModel(cor_dat, discrete_model, rate.cat, dual, collapse)
    index.disc[index.disc == 0] <- NA
  }else{
    index.disc <- discrete_model
    index.disc[index.disc == 0] <- NA
  }
  if(class(continuous_model)[1] == "character"){
    index.cont <- getOUParamStructure(continuous_model, nStates, rate.cat, null.model)
  }else{
    continuous_model[continuous_model == 0] <- NA
    index.cont <- continuous_model
  }
  if(dim(index.disc)[2] > dim(index.cont)[2]){
    stop("Not all of your discrete states have OU parameters associated with them. Please check that your discrete index matrix matches your continuous index matrix.")
  }
  if(dim(index.cont)[2] > dim(index.disc)[2]){
    stop("You have specified more OU parameters than there are states in the discrete process. Please check that your discrete index matrix matches your continuous index matrix.")
  }
  if(class(root.p[1]) != "character"){
    if(dim(index.disc)[2] != length(root.p)){
      stop("You have entered a custom root prior whose length does not equal the number of states in your discrete model.")
    }
  }
  
  if(is.null(lb_continuous_model)){
    # the lower limit of alpha is defined as a halflife of 10000% of the max tree height
    # the lower limit of sigma is defined 10 times less than alpha
    # the lower limit of optim is defined 10 times lower than the minimum observation
    if(any(is.na(lb_continuous_model[1,]))){
      lb.alpha = 1e-10
    }else{
      lb.alpha = 1e-10
    }
    lb.sigma = 1e-10
    lb.optim = min(data[, 1+nCol+1])/10 
    lb_continuous_model=c(lb.alpha,lb.sigma,lb.optim)
  }
  if(is.null(ub_continuous_model)){
    # the upper limit of alpha is defined as a halflife of 1% of the max tree height
    # the upper limit of sigma is defined 10 times more than alpha
    # the upper limit of optim is defined 10 times more than the maximum observation
    ub.alpha = log(2)/(0.01 * Tmax)
    ub.sigma = ub.alpha
    ub.optim = max(data[, 1+nCol+1])*10 
    ub_continuous_model=c(ub.alpha,ub.sigma,ub.optim)
  }
  if(is.null(lb_discrete_model)){
    # the minimum dwell time is defined as 10000 times the max tree height
    lb_discrete_model = 1/(Tmax*10000)
  }
  if(is.null(ub_discrete_model)){
    ub_discrete_model = 1/(Tmax*0.0001)
  }
  #Ensures that weird root state probabilities that do not sum to 1 are input:
  if(!is.null(root.p)){
    if(!is.character(root.p)){
      root.p <- root.p/sum(root.p)
    }
  }
  time_slice <- Tmax+1
  # the number of parameters for each process
  n_p_trans <- max(index.disc, na.rm = TRUE)
  n_p_alpha <- length(unique(na.omit(index.cont[1,])))
  n_p_sigma <- length(unique(na.omit(index.cont[2,])))
  n_p_theta <- length(unique(na.omit(index.cont[3,])))
  n_p <- n_p_trans + n_p_alpha + n_p_sigma + n_p_theta
  
  # an internal data structure (internodes liks matrix) for the dev function
  edge_liks_list <- getEdgeLiks(phy, hOUwie.dat$data.cor, nStates, rate.cat, time_slice)
  
  # default MLE search options
  if(is.null(opts)){
    if(optimizer == "nlopt_ln"){
      opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.25)
    }
    if(optimizer == "nlopt_gn"){
      opts <- list("algorithm"="NLOPT_GN_DIRECT_L", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.25)
    }
    if(optimizer == "sann"){
      opts <- list(max.call=1000, smooth=FALSE)
    }
  }
  # a global matrix to contain likelihoods so that identical parameters return identical likelihoods
  if(is.null(opts$maxeval) | is.null(opts$max.call)){
    max.its <- 1000
  }else{
    max.its <- as.numeric(opts$maxeval)
  }
  setDTthreads(threads=1)
  tmp.df <- data.frame(matrix(c(0, rep(1e5, n_p)), byrow = TRUE, ncol = n_p+1, nrow = max.its))
  global_liks_mat <- as.data.table(tmp.df)
  
  # p is organized into 2 groups with the first set being corHMM and the second set being OUwie
  # organized as c(trans.rt, alpha, sigma.sq, theta)
  # evaluate likelihood
  if(!is.null(p)){
    if(negative_values){
      p[(n_p_trans + n_p_alpha + n_p_sigma + 1):n_p] <- p[(n_p_trans + n_p_alpha + n_p_sigma + 1):n_p] + 50 
    }
    if(!quiet){
      cat("Calculating likelihood from a set of fixed parameters.\n")
      print(p)
    }
    if(max(index.cont, na.rm = TRUE) + max(index.disc, na.rm = TRUE) != length(p)){
      message <- paste0("The number of parameters does not match the number required by the model structure. You have supplied ", length(p), ", but the model structure requires ", max(index.cont, na.rm = TRUE) + max(index.disc, na.rm = TRUE), ".")
      stop(message, call. = FALSE)
    }
    out<-NULL
    pars <- out$solution <- log(p)
    # out$objective <- hOUwie.dev(p = log(p), phy=phy, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p,edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths, sample_tips=sample_tips, sample_nodes=sample_nodes, split.liks=FALSE)
  }else{
    out<-NULL
    lower = log(c(rep(lb_discrete_model, n_p_trans), 
                  rep(lb_continuous_model[1], length(unique(na.omit(index.cont[1,])))), 
                  rep(lb_continuous_model[2], length(unique(na.omit(index.cont[2,])))), 
                  rep(lb_continuous_model[3], length(unique(na.omit(index.cont[3,]))))))
    upper = log(c(rep(ub_discrete_model, n_p_trans), 
                  rep(ub_continuous_model[1], length(unique(na.omit(index.cont[1,])))), 
                  rep(ub_continuous_model[2], length(unique(na.omit(index.cont[2,])))), 
                  rep(ub_continuous_model[3], length(unique(na.omit(index.cont[3,]))))))
    # cat(c("TotalLnLik", "DiscLnLik", "ContLnLik"), "\n")
    # check for user input initial parameters 
    if(is.null(ip)){
      if(rate.cat > 1){
        bin_index <- cut(hOUwie.dat$data.ou[,3], rate.cat, labels = FALSE)
        combos <- expand.grid(1:max(hOUwie.dat$data.cor[,2]), 1:rate.cat)
        disc_tips <- vector("numeric", length(phy$tip.label))
        for(i in 1:dim(combos)[1]){
          disc_tips[hOUwie.dat$data.cor[,2] == combos[i,1] & bin_index == combos[i,2]] <- i
        }
      }else{
        disc_tips <- hOUwie.dat$data.cor[,2]
      }
      starts.alpha <- rep(log(2)/Tmax, n_p_alpha)
      # starts.sigma <- rep(var(hOUwie.dat$data.ou[,3]), n_p_sigma)
      starts.sigma <- rep(log(2)/Tmax, n_p_sigma)
      start.theta <- getIP.theta(hOUwie.dat$data.ou[,3], disc_tips, index.cont[3,])
      start.cor <- rep(10/sum(phy$edge.length), n_p_trans)
      starts.basic = c(start.cor, starts.alpha, starts.sigma, start.theta)
      starts <- starts.basic
    }else{
      if(negative_values){
        ip[(n_p_trans + n_p_alpha + n_p_sigma + 1):n_p] <- ip[(n_p_trans + n_p_alpha + n_p_sigma + 1):n_p] + 50 
      }
      starts <- ip
    }
    if(!quiet){
      cat("Starting a thorough search with", nSim, "simmaps using the", optimizer, "optimization protocol...\n")
    }
    multiple_starts <- generateMultiStarting(starts, index.disc, index.cont, n_starts, exp(lower), exp(upper))
    if(length(grep("nlopt", optimizer)) == 1){
      # out = nloptr(x0=log(starts), eval_f=hOUwie.dev, lb=lower, ub=upper, opts=opts,
      #              phy=phy, data=hOUwie.dat$data.ou,
      #              rate.cat=rate.cat, mserr=mserr,
      #              index.disc=index.disc, index.cont=index.cont, root.p=root.p,
      #              edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths,
      #              sample_tips=sample_tips, split.liks=FALSE)
      multi_out <- mclapply(multiple_starts, function(x) nloptr(x0=log(x), eval_f=hOUwie.dev, lb=lower, ub=upper, opts=opts, phy=phy, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr,index.disc=index.disc, index.cont=index.cont, root.p=root.p,edge_liks_list=edge_liks_list, nSim=nSim, all.paths=all.paths, sample_tips=sample_tips, sample_nodes=sample_nodes, adaptive_sampling=adaptive_sampling, split.liks=FALSE, global_liks_mat=global_liks_mat, diagn_msg=diagn_msg), mc.cores = ncores)
      multi_logliks <- unlist(lapply(multi_out, function(x) x$objective))
      if(any(-multi_logliks > 1e10) | any(is.null(multi_logliks))){
        cat("\nIt appears that an optimization failed. Removing failed optimizations from final output.\n")
        failed_optimizations <- which(-multi_logliks > 1e10 | is.null(multi_logliks))
        multi_logliks <- multi_logliks[-failed_optimizations]
        multi_out <- multi_out[-failed_optimizations]
        if(length(multi_out) == 0){
          return(NULL)
        }
      }
      out <- multi_out[[which.min(multi_logliks)]]
      pars <- out$solution
    }
    if(length(grep("sann", optimizer)) == 1){
      # out = GenSA(par=log(starts), fn=hOUwie.dev, lower=lower, upper=upper, control=opts, 
      #              phy=phy, data=hOUwie.dat$data.ou, 
      #              rate.cat=rate.cat, mserr=mserr, 
      #              index.disc=index.disc, index.cont=index.cont, root.p=root.p,
      #              edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths, 
      #              sample_tips=sample_tips, split.liks=FALSE)
      multi_out <- mclapply(multiple_starts, function(x) GenSA(par=log(x), fn=hOUwie.dev, lower=lower, upper=upper, control=opts, phy=phy, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, edge_liks_list=edge_liks_list, nSim=nSim, all.paths=all.paths, sample_tips=sample_tips, sample_nodes=sample_nodes, adaptive_sampling=adaptive_sampling, split.liks=FALSE), global_liks_mat=global_liks_mat, diagn_msg=diagn_msg, mc.cores = ncores)
      multi_logliks <- unlist(lapply(multi_out, function(x) x$value))
      out <- multi_out[[which.min(multi_logliks)]]
      pars <- out$par
    }
  }
  # preparing output
  liks_houwie <- hOUwie.dev(p = pars, phy=phy, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, edge_liks_list=edge_liks_list, nSim=nSim, all.paths=all.paths, sample_tips=sample_tips, sample_nodes=sample_nodes, adaptive_sampling=adaptive_sampling, split.liks=TRUE, global_liks_mat=global_liks_mat, diagn_msg=FALSE)
  houwie_obj <- getHouwieObj(liks_houwie, pars=exp(pars), phy=phy, data=data, hOUwie.dat=hOUwie.dat, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, nSim=nSim, sample_tips=sample_tips, sample_nodes=sample_nodes, adaptive_sampling=adaptive_sampling, nStates=nStates, discrete_model=discrete_model, continuous_model=continuous_model, time_slice=time_slice, root.station=root.station, get.root.theta=get.root.theta,lb_discrete_model,ub_discrete_model,lb_continuous_model,ub_continuous_model, ip=ip, opts=opts, quiet=quiet, negative_values=negative_values)
  # adding independent model if included
  # if(is.null(p)){
  #   liks_indep <- hOUwie.dev(p = log(starts), phy=phy, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths, sample_tips=sample_tips, split.liks=TRUE)
  #   houwie_obj$init_model <- getHouwieObj(liks_indep, pars=starts, phy=phy, data=data, hOUwie.dat=hOUwie.dat, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, nSim=nSim, sample_tips=sample_tips, nStates=nStates, discrete_model=discrete_model, continuous_model=continuous_model, time_slice=time_slice, root.station=root.station, get.root.theta=get.root.theta,lb_discrete_model,ub_discrete_model,lb_continuous_model,ub_continuous_model, ip=ip, opts=opts, quiet=quiet)
  # }
  # conducting ancestra state resconstruction
  if(recon){
    houwie_recon <- hOUwie.recon(houwie_obj, nodes)
    houwie_obj$recon <- houwie_recon
  }
  houwie_obj$all_disc_liks <- liks_houwie$llik_discrete
  houwie_obj$all_cont_liks <- liks_houwie$llik_continuous
  houwie_obj$simmaps <- lapply(liks_houwie$simmaps, correct_map_edges)
  houwie_obj$global_liks_mat <- global_liks_mat
  end_time <- Sys.time()
  run_time <- end_time - start_time
  houwie_obj$run_time <- run_time
  units(houwie_obj$run_time) <- "mins"
  return(houwie_obj)
}

hOUwie.fixed <- function(simmaps, data, rate.cat, discrete_model, continuous_model, null.model=FALSE, root.p="yang", dual = FALSE, collapse = TRUE, root.station=FALSE, get.root.theta=FALSE, mserr = "none", lb_discrete_model=NULL, ub_discrete_model=NULL, lb_continuous_model=NULL, ub_continuous_model=NULL, recon=FALSE, nodes="internal", p=NULL, ip=NULL, optimizer="nlopt_ln", opts=NULL, quiet=FALSE, sample_tips=FALSE, sample_nodes=TRUE, adaptive_sampling=FALSE, diagn_msg=FALSE, make_numeric = TRUE, n_starts = 1, ncores = 1){
  start_time <- Sys.time()
  # if the data has negative values, shift it right - we will shift it back later
  negative_values <- FALSE
  if(mserr == "none"){
    if(any(data[,dim(data)[2]] < 0)){
      cat("Negative values detected... adding 50 to the trait mean for optimization purposes\n")
      negative_values <- TRUE
      data[,dim(data)[2]] <- data[,dim(data)[2]] + 50
    }
  }else{
    if(any(data[,dim(data)[2]-1] < 0)){
      cat("Negative values detected... adding 50 to the trait mean for optimization purposes\n")
      negative_values <- TRUE
      data[,dim(data)[2]-1] <- data[,dim(data)[2]-1] + 50
    }
  }
  if(ncores > n_starts){
    cat("You have specified more cores are to be used than the number of starts. Setting ncores to be equal to the number of optimizations.\n")
    ncores <- n_starts
  }
  
  # organize the data
  hOUwie.dat <- organizeHOUwieDat(data, mserr, collapse)
  observed_traits <- hOUwie.dat$ObservedTraits
  names(observed_traits) <- 1:length(observed_traits)
  cat("\nUsing the following legend:\n")
  print(observed_traits)
  if(make_numeric){
    simmaps <- lapply(simmaps, function(x) makeMapEdgesNumeric(x, observed_traits))
  }
  nStates <- as.numeric(max(hOUwie.dat$data.cor[,2]))
  nCol <- dim(data)[2] - ifelse(mserr == "none", 2, 3)
  Tmax <- max(branching.times(simmaps[[1]]))
  all.paths <- lapply(1:(Nnode(simmaps[[1]]) + Ntip(simmaps[[1]])), function(x) getPathToRoot(simmaps[[1]], x))
  
  if(class(discrete_model)[1] == "character"){
    index.disc <- getDiscreteModel(hOUwie.dat$data.cor, discrete_model, rate.cat, dual, collapse)
    index.disc[index.disc == 0] <- NA
  }else{
    index.disc <- discrete_model
    index.disc[index.disc == 0] <- NA
  }
  if(class(continuous_model)[1] == "character"){
    index.cont <- getOUParamStructure(continuous_model, nStates, rate.cat, null.model)
  }else{
    continuous_model[continuous_model == 0] <- NA
    index.cont <- continuous_model
  }
  if(dim(index.disc)[2] > dim(index.cont)[2]){
    stop("Not all of your discrete states have OU parameters associated with them. Please check that your discrete index matrix matches your continuous index matrix.")
  }
  if(dim(index.cont)[2] > dim(index.disc)[2]){
    stop("You have specified more OU parameters than there are states in the discrete process. Please check that your discrete index matrix matches your continuous index matrix.")
  }
  if(class(root.p[1]) != "character"){
    if(dim(index.disc)[2] != length(root.p)){
      stop("You have entered a custom root prior whose length does not equal the number of states in your discrete model.")
    }
  }
  
  if(is.null(lb_continuous_model)){
    # the lower limit of alpha is defined as a halflife of 10000% of the max tree height
    # the lower limit of sigma is defined 10 times less than alpha
    # the lower limit of optim is defined 10 times lower than the minimum observation
    if(any(is.na(lb_continuous_model[1,]))){
      lb.alpha = 1e-10
    }else{
      lb.alpha = 1e-10
    }
    lb.sigma = 1e-10
    lb.optim = min(data[, 1+nCol+1])/10 
    lb_continuous_model=c(lb.alpha,lb.sigma,lb.optim)
  }
  if(is.null(ub_continuous_model)){
    # the upper limit of alpha is defined as a halflife of 1% of the max tree height
    # the upper limit of sigma is defined 10 times more than alpha
    # the upper limit of optim is defined 10 times more than the maximum observation
    ub.alpha = log(2)/(0.01 * Tmax)
    ub.sigma = ub.alpha
    ub.optim = max(data[, 1+nCol+1])*10 
    ub_continuous_model=c(ub.alpha,ub.sigma,ub.optim)
  }
  if(is.null(lb_discrete_model)){
    # the minimum dwell time is defined as 100 times the max tree height
    lb_discrete_model = 1/(Tmax*100)
  }
  if(is.null(ub_discrete_model)){
    ub_discrete_model = 1/(Tmax*0.01)
  }
  #Ensures that weird root state probabilities that do not sum to 1 are input:
  if(!is.null(root.p)){
    if(!is.character(root.p)){
      root.p <- root.p/sum(root.p)
    }
  }
  time_slice <- Tmax+1
  # the number of parameters for each process
  n_p_trans <- max(index.disc, na.rm = TRUE)
  n_p_alpha <- length(unique(na.omit(index.cont[1,])))
  n_p_sigma <- length(unique(na.omit(index.cont[2,])))
  n_p_theta <- length(unique(na.omit(index.cont[3,])))
  n_p <- n_p_trans + n_p_alpha + n_p_sigma + n_p_theta
  
  # an internal data structure (internodes liks matrix) for the dev function
  edge_liks_list <- getEdgeLiks(simmaps[[1]], hOUwie.dat$data.cor, nStates, rate.cat, time_slice)
  
  # default MLE search options
  if(is.null(opts)){
    if(optimizer == "nlopt_ln"){
      opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.25)
    }
    if(optimizer == "nlopt_gn"){
      opts <- list("algorithm"="NLOPT_GN_DIRECT_L", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.25)
    }
    if(optimizer == "sann"){
      opts <- list(max.call=1000, smooth=FALSE)
    }
  }
  # a global matrix to contain likelihoods so that identical parameters return identical likelihoods
  if(is.null(opts$maxeval) | is.null(opts$max.call)){
    max.its <- 1000
  }else{
    max.its <- as.numeric(opts$maxeval)
  }
  setDTthreads(threads=1)
  tmp.df <- data.frame(matrix(c(0, rep(1e5, n_p)), byrow = TRUE, ncol = n_p+1, nrow = max.its))
  global_liks_mat <- as.data.table(tmp.df)
  
  # p is organized into 2 groups with the first set being corHMM and the second set being OUwie
  # organized as c(trans.rt, alpha, sigma.sq, theta)
  # evaluate likelihood
  if(!is.null(p)){
    if(negative_values){
      p[(n_p_trans + n_p_alpha + n_p_sigma + 1):n_p] <- p[(n_p_trans + n_p_alpha + n_p_sigma + 1):n_p] + 50 
    }
    if(!quiet){
      cat("Calculating likelihood from a set of fixed parameters.\n")
      print(p)
    }
    if(max(index.cont, na.rm = TRUE) + max(index.disc, na.rm = TRUE) != length(p)){
      message <- paste0("The number of parameters does not match the number required by the model structure. You have supplied ", length(p), ", but the model structure requires ", max(index.cont, na.rm = TRUE) + max(index.disc, na.rm = TRUE), ".")
      stop(message, call. = FALSE)
    }
    out<-NULL
    pars <- out$solution <- log(p)
    # out$objective <- hOUwie.dev(p = log(p), phy=phy, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p,edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths, sample_tips=sample_tips, sample_nodes=sample_nodes, split.liks=FALSE)
  }else{
    out<-NULL
    lower = log(c(rep(lb_discrete_model, n_p_trans), 
                  rep(lb_continuous_model[1], length(unique(na.omit(index.cont[1,])))), 
                  rep(lb_continuous_model[2], length(unique(na.omit(index.cont[2,])))), 
                  rep(lb_continuous_model[3], length(unique(na.omit(index.cont[3,]))))))
    upper = log(c(rep(ub_discrete_model, n_p_trans), 
                  rep(ub_continuous_model[1], length(unique(na.omit(index.cont[1,])))), 
                  rep(ub_continuous_model[2], length(unique(na.omit(index.cont[2,])))), 
                  rep(ub_continuous_model[3], length(unique(na.omit(index.cont[3,]))))))
    # cat(c("TotalLnLik", "DiscLnLik", "ContLnLik"), "\n")
    # check for user input initial parameters 
    if(is.null(ip)){
      if(rate.cat > 1){
        bin_index <- cut(hOUwie.dat$data.ou[,3], rate.cat, labels = FALSE)
        combos <- expand.grid(1:max(hOUwie.dat$data.cor[,2]), 1:rate.cat)
        disc_tips <- vector("numeric", length(simmaps[[1]]$tip.label))
        for(i in 1:dim(combos)[1]){
          disc_tips[hOUwie.dat$data.cor[,2] == combos[i,1] & bin_index == combos[i,2]] <- i
        }
      }else{
        disc_tips <- hOUwie.dat$data.cor[,2]
      }
      starts.alpha <- rep(log(2)/Tmax, n_p_alpha)
      # starts.sigma <- rep(var(hOUwie.dat$data.ou[,3]), n_p_sigma)
      starts.sigma <- rep(log(2)/Tmax, n_p_sigma)
      start.theta <- getIP.theta(hOUwie.dat$data.ou[,3], disc_tips, index.cont[3,])
      start.cor <- rep(10/sum(simmaps[[1]]$edge.length), n_p_trans)
      starts.basic = c(start.cor, starts.alpha, starts.sigma, start.theta)
      starts <- starts.basic
    }else{
      if(negative_values){
        ip[(n_p_trans + n_p_alpha + n_p_sigma + 1):n_p] <- ip[(n_p_trans + n_p_alpha + n_p_sigma + 1):n_p] + 50 
      }
      starts <- ip
    }
    if(!quiet){
      cat("Starting a thorough search using the", optimizer, "optimization protocol...\n")
    }
    multiple_starts <- generateMultiStarting(starts, index.disc, index.cont, n_starts, exp(lower), exp(upper))
    if(length(grep("nlopt", optimizer)) == 1){
      # out = nloptr(x0=log(starts), eval_f=hOUwie.dev, lb=lower, ub=upper, opts=opts,
      #              phy=phy, data=hOUwie.dat$data.ou,
      #              rate.cat=rate.cat, mserr=mserr,
      #              index.disc=index.disc, index.cont=index.cont, root.p=root.p,
      #              edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths,
      #              sample_tips=sample_tips, split.liks=FALSE)
      multi_out <- mclapply(multiple_starts, function(x) nloptr(x0=log(x), eval_f=hOUwie.fixed.dev, lb=lower, ub=upper, opts=opts, simmaps=simmaps, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr,index.disc=index.disc, index.cont=index.cont, root.p=root.p,edge_liks_list=edge_liks_list, all.paths=all.paths, sample_tips=sample_tips, sample_nodes=sample_nodes, adaptive_sampling=adaptive_sampling, split.liks=FALSE, global_liks_mat=global_liks_mat, diagn_msg=diagn_msg), mc.cores = ncores)
      multi_logliks <- unlist(lapply(multi_out, function(x) x$objective))
      search_summary <- c(best_loglik = -min(multi_logliks), mean_loglik = -log(mean(exp(multi_logliks))), sd_logliks = log(sd(exp(multi_logliks))))
      if(!quiet){
        cat("\nOptimization complete. Optimization summary:\n")
        print(search_summary)
      }
      out <- multi_out[[which.min(multi_logliks)]]
      pars <- out$solution
    }
    if(length(grep("sann", optimizer)) == 1){
      # out = GenSA(par=log(starts), fn=hOUwie.dev, lower=lower, upper=upper, control=opts, 
      #              phy=phy, data=hOUwie.dat$data.ou, 
      #              rate.cat=rate.cat, mserr=mserr, 
      #              index.disc=index.disc, index.cont=index.cont, root.p=root.p,
      #              edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths, 
      #              sample_tips=sample_tips, split.liks=FALSE)
      multi_out <- mclapply(multiple_starts, function(x) GenSA(par=log(x), fn=hOUwie.fixed.dev, lower=lower, upper=upper, control=opts, simmaps=simmaps, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, edge_liks_list=edge_liks_list, all.paths=all.paths, sample_tips=sample_tips, sample_nodes=sample_nodes, adaptive_sampling=adaptive_sampling, split.liks=FALSE, diagn_msg=diagn_msg), global_liks_mat=global_liks_mat, mc.cores = ncores)
      multi_logliks <- unlist(lapply(multi_out, function(x) x$value))
      search_summary <- c(best_loglik = -min(multi_logliks), mean_loglik = -mean(multi_logliks), sd_logliks = sd(multi_logliks))
      if(!quiet){
        cat("Optimization complete. Optimization summary:")
        print(search_summary)
      }
      out <- multi_out[[which.min(multi_logliks)]]
      pars <- out$par
    }
  }
  # preparing output
  liks_houwie <- hOUwie.fixed.dev(p = pars, simmaps=simmaps, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, edge_liks_list=edge_liks_list, all.paths=all.paths, sample_tips=sample_tips, sample_nodes=sample_nodes, adaptive_sampling=adaptive_sampling, split.liks=TRUE, global_liks_mat=global_liks_mat)
  houwie_obj <- getHouwieObj(liks_houwie, pars=exp(pars), phy=simmaps[[1]], data=data, hOUwie.dat=hOUwie.dat, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, nSim=NULL, sample_tips=sample_tips, sample_nodes=sample_nodes, adaptive_sampling=adaptive_sampling, nStates=nStates, discrete_model=discrete_model, continuous_model=continuous_model, time_slice=time_slice, root.station=root.station, get.root.theta=get.root.theta,lb_discrete_model,ub_discrete_model,lb_continuous_model,ub_continuous_model, ip=ip, opts=opts, quiet=quiet, negative_values=negative_values)
  # adding independent model if included
  # if(is.null(p)){
  #   liks_indep <- hOUwie.dev(p = log(starts), phy=phy, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, edge_liks_list=edge_liks_list, nSim=nSim, tip.paths=tip.paths, sample_tips=sample_tips, split.liks=TRUE)
  #   houwie_obj$init_model <- getHouwieObj(liks_indep, pars=starts, phy=phy, data=data, hOUwie.dat=hOUwie.dat, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, nSim=nSim, sample_tips=sample_tips, nStates=nStates, discrete_model=discrete_model, continuous_model=continuous_model, time_slice=time_slice, root.station=root.station, get.root.theta=get.root.theta,lb_discrete_model,ub_discrete_model,lb_continuous_model,ub_continuous_model, ip=ip, opts=opts, quiet=quiet)
  # }
  # conducting ancestra state resconstruction
  if(recon){
    houwie_recon <- hOUwie.recon(houwie_obj, nodes)
    houwie_obj$recon <- houwie_recon
  }
  houwie_obj$all_disc_liks <- liks_houwie$llik_discrete
  houwie_obj$all_cont_liks <- liks_houwie$llik_continuous
  houwie_obj$simmaps <- lapply(liks_houwie$simmaps, correct_map_edges)   
  houwie_obj$global_liks_mat <- global_liks_mat
  end_time <- Sys.time()
  run_time <- end_time - start_time
  houwie_obj$run_time <- run_time
  units(houwie_obj$run_time) <- "mins"
  return(houwie_obj)
}


hOUwie.recon <- function(houwie_obj, nodes="all"){
  # if the class is houwie_obj
  phy <- houwie_obj$phy
  hOUwie.dat <- houwie_obj$hOUwie.dat
  root.p <- houwie_obj$root.p
  mserr <- houwie_obj$mserr
  rate.cat <- houwie_obj$rate.cat
  index.cont <- houwie_obj$index.cont
  index.disc <- houwie_obj$index.disc
  p <- houwie_obj$p
  time_slice <- houwie_obj$time_slice
  state_names <- colnames(houwie_obj$solution.disc)
  sample_tips <- houwie_obj$sample_tips
  sample_nodes <- houwie_obj$sample_nodes
  adaptive_sampling <- houwie_obj$adaptive_sampling
  nSim <- houwie_obj$nSim
  # organize the data
  phy <- reorder.phylo(phy, "pruningwise")
  nTip <- length(phy$tip.label)
  Tmax <- max(branching.times(phy))
  nStates <- as.numeric(max(hOUwie.dat$data.cor[,2]))
  all.paths <- lapply(1:(Nnode(phy) + Ntip(phy)), function(x) getPathToRoot(phy, x))
  # an internal data structure (internodes liks matrix) for the dev function
  edge_liks_list <- getEdgeLiks(phy, hOUwie.dat$data.cor, nStates, rate.cat, time_slice)
  if(is.character(nodes[1])){
    if(nodes == "internal"){
      nodes_to_fix <- min(phy$edge[,1]):max(phy$edge[,1])
    }
    if(nodes == "external"){
      nodes_to_fix <- 1:nTip
    }
    if(nodes  == "all"){
      nodes_to_fix <- 1:max(phy$edge[,1])
    }
  }else{
    nodes_to_fix <- nodes
  }
  recon_matrix <- matrix(-Inf, length(nodes_to_fix), nStates * rate.cat, dimnames = list(nodes_to_fix, state_names))
  for(i in 1:length(nodes_to_fix)){
    cat("\rState reconstruction for node", i, "of", length(nodes_to_fix), "...")
    node_i <- nodes_to_fix[i]
    anc_edges_to_fix <- which(phy$edge[,1] == node_i)
    dec_edges_to_fix <- which(phy$edge[,2] == node_i)
    if(node_i <= nTip){
      possible_states <- which(edge_liks_list[[dec_edges_to_fix]][1,] == 1)
    }else{
      possible_states <- 1:(nStates*rate.cat)
    }
    # check_unreasonable_recon <- TRUE
    # while(check_unreasonable_recon){
    for(state_j in 1:(nStates*rate.cat)){
      if(!state_j %in% possible_states){
        next
      }
      edge_liks_list_i <- edge_liks_list
      fix_vector <- numeric(nStates * rate.cat)
      fix_vector[state_j] <- 1
      for(k in dec_edges_to_fix){
        edge_liks_list_i[[k]][1,] <- fix_vector
      }
      for(k in anc_edges_to_fix){
        last_row <- dim(edge_liks_list_i[[k]])[1]
        edge_liks_list_i[[k]][last_row,] <- fix_vector
      }
      fixed_loglik <- -hOUwie.dev(p = log(p), phy=phy, data=hOUwie.dat$data.ou, rate.cat=rate.cat, mserr=mserr, index.disc=index.disc, index.cont=index.cont, root.p=root.p, edge_liks_list=edge_liks_list_i, nSim=nSim, all.paths=all.paths, sample_tips=sample_tips, sample_nodes=sample_nodes, adaptive_sampling=adaptive_sampling, split.liks=FALSE, global_liks_mat=NULL)
      recon_matrix[i, state_j] <- fixed_loglik
    }
    recon_loglik <- max(recon_matrix[i, ]) + log(sum(exp(recon_matrix[i, ] - max(recon_matrix[i, ]))))
    # check_unreasonable_recon <- recon_loglik < houwie_obj$loglik * 1.05
    # }
  }
  recon_matrix <- t(apply(recon_matrix, 1, function(x) x - max(x)))
  recon_matrix <- round(exp(recon_matrix)/rowSums(exp(recon_matrix)), 10)
  return(recon_matrix)
}

# simulate a hOUwie model
hOUwie.sim <- function(phy, Q, root.freqs, alpha, sigma.sq, theta0, theta){
  # simulate an Mk dataset
  dat.cor <- simCharacterHistory(phy, Q, root.freqs)
  while(!all(1:dim(Q)[1] %in% dat.cor$TipStates)){
    dat.cor <- simCharacterHistory(phy, Q, root.freqs)
  }
  map <- getMapFromNode(phy, dat.cor$TipStates, dat.cor$NodeStates, 0.5)
  simmap <- getMapFromSubstHistory(list(map), phy)[[1]]
  # dat.cor <- data.frame(sp=names(dat.cor), d=dat.cor)
  # simulate a stochastic map with true Q
  # simmap <- corHMM:::makeSimmap(phy, dat.cor, Q, 1, nSim=nMap)
  # lik <- corHMM:::getSimmapLik(simmap, Q)
  # simulate the ou dataset
  dat.ou <- OUwie.sim(simmap, simmap.tree = TRUE, alpha = alpha, sigma.sq = sigma.sq, theta0 = theta0, theta = theta)
  # dat.ou <- OUwie.sim(simmap, simmap.tree = TRUE, alpha = alpha, sigma.sq = sig2, theta0 = theta0, theta = theta)
  # return true params and data
  data <- data.frame(sp = phy$tip.label, reg = dat.cor$TipStates, x = dat.ou$X)
  rownames(data) <- NULL
  return(list(data = data, simmap = simmap))
}

# rerun a set of completed models with the best current maps
hOUwie.thorough <- function(model.list, ncores=1){
  # seprate models by their rate class
  rate_cat_vector <- unlist(lapply(model.list, "[[", "rate.cat"))
  nSim_vector <- unlist(lapply(model.list, "[[", "nSim"))
  if(length(unique(nSim_vector)) > 1){
    stop("You have different amounts of simmaps for each model estimation. These should all be the same.")
  }else{
    nSim <- unique(nSim_vector)
  }
  all_rate_cats <- unique(rate_cat_vector)
  model_set_separated <- sapply(all_rate_cats, function(x) model.list[x == rate_cat_vector])
  # for each rate class run hOUwie summarize the maps
  model_avg_pars <- lapply(model_set_separated, getModelAvgParams)
  new_res <- list()
  for(i in 1:length(model_set_separated)){
    cat("Preparing rate category", i, "for hOUwie thorough...\n")
    all_maps <- do.call(c, lapply(model_set_separated[[i]], "[[", "simmaps"))
    map_id_list <- unlist(lapply(all_maps, function(x) paste0(names(unlist(x$maps)), collapse = "")))
    map_lik_list <- unlist(lapply(model_set_separated[[i]], function(x) x$all_cont_liks + x$all_disc_liks))
    unique_map_index <- !duplicated(map_id_list)
    unique_map_lik_list <- map_lik_list[unique_map_index]
    uniqie_all_maps <- all_maps[unique_map_index]
    sorted_index <- sort(unique_map_lik_list, decreasing = TRUE, index.return = TRUE)$ix
    new_maps <- uniqie_all_maps[sorted_index[1:nSim]]
    new_res[[i]] <- mclapply(model_set_separated[[i]], function(x) runSingleThorough(x, new_maps, model_avg_pars[[i]]), mc.cores = ncores)
  }
  out <- do.call(c, new_res)
  names(out) <- names(model.list)
  return(out)
}

##### Utility exported functions ##### 
hOUwie.walk <- function(houwie_obj, delta=2, nsteps=1000, print_freq=50, lower_bound=0, upper_bound=Inf, adjust_width_interval=100, badval=1e9, sd_vector=NULL, debug=FALSE, restart_after=50){
  
  phy <- houwie_obj$phy
  root.p <- houwie_obj$root.p
  mserr <- houwie_obj$mserr
  rate.cat <- houwie_obj$rate.cat
  index.cont <- houwie_obj$index.cont
  index.disc <- houwie_obj$index.disc
  time_slice <- houwie_obj$time_slice
  state_names <- colnames(houwie_obj$solution.disc)
  sample_tips <- houwie_obj$sample_tips
  sample_nodes <- houwie_obj$sample_nodes
  adaptive_sampling <- houwie_obj$adaptive_sampling
  nSim <- houwie_obj$nSim
  data <- houwie_obj$data
  # set up dentist
  best_par <- (houwie_obj$p)
  n_p_trans <- max(index.disc, na.rm = TRUE)
  n_p_alpha <- length(unique(na.omit(index.cont[1,])))
  n_p_sigma <- length(unique(na.omit(index.cont[2,])))
  n_p_theta <- length(unique(na.omit(index.cont[3,])))
  names(best_par) <- c(paste0("rate", "_", 1:n_p_trans), paste0("alpha", "_", 1:n_p_alpha), paste0("sigma2", "_", 1:n_p_sigma), paste0("theta", "_", 1:n_p_theta))
  best_neglnL <- -houwie_obj$loglik
  
  houwie_to_run <- function(par){
    fixed_loglik <- hOUwie(p = par, phy=phy, data=houwie_obj$data, rate.cat=rate.cat, mserr=mserr, discrete_model = index.disc, continuous_model = index.cont, root.p=root.p, nSim=nSim, sample_tips=sample_tips, sample_nodes=sample_nodes, adaptive_sampling=adaptive_sampling, quiet = TRUE)$loglik
    return(-fixed_loglik)
  }
  
  dented_results <- dent_walk(par=best_par, fn=houwie_to_run, best_neglnL=best_neglnL, delta=delta, nsteps=nsteps, print_freq=print_freq, lower_bound=lower_bound, upper_bound=upper_bound, adjust_width_interval=adjust_width_interval, badval=badval, sd_vector=sd_vector, debug=debug, restart_after=restart_after)
  
  return(dented_results)
}

getModelTable <- function(model.list, type="AIC"){
  # checks
  if(class(model.list) != "list"){
    stop("Input object must be of class list with each element as a separet fit model to the same dataset.", call. = FALSE)
  }
  if(!all(unlist(lapply(model.list, function(x) class(x))) == "houwie")){
    warning("Not all models are of class houwie. These have been removed.")
    model.list <- model.list[unlist(lapply(model.list, function(x) class(x)) == "houwie")]
  }
  if(var(unlist(lapply(model.list, function(x) dim(x$data)[1]))) != 0){
    stop("The number of rows in your data are not the same for all models. Models should not be compared if they are not evaluating the same dataset.", call.=FALSE)
  }
  if(length(model.list) == 1){
    stop("Two or models are needed to conduct model averaging.", call. = FALSE)
  }
  
  ParCount <- unlist(lapply(model.list, function(x) x$param.count))
  nTip <- length(model.list[[1]]$phy$tip.label)
  AIC <- simplify2array(lapply(model.list, "[[", type))
  dAIC <- AIC - min(AIC)
  AICwt <- exp(-0.5 * dAIC)/sum(exp(-0.5 * dAIC))
  LogLik <- simplify2array(lapply(model.list, "[[", "loglik"))
  DiscLik <- simplify2array(lapply(model.list, "[[", "DiscLik"))
  ContLik <- simplify2array(lapply(model.list, "[[", "ContLik"))
  model_table <- data.frame(np = ParCount, lnLik = LogLik, DiscLik=DiscLik, ContLik=ContLik, AIC = AIC, dAIC = dAIC, AICwt = AICwt)
  colnames(model_table) <- gsub("AIC", type, colnames(model_table))
  return(model_table)
}

getModelAvgParams <- function(model.list, BM_alpha_treatment="zero", force=TRUE){
  if(any(unlist(lapply(model.list, class))!="houwie")){
    warning("Some of the input models are not of class houwie, these have been removed.")
    model.list <- model.list[which(unlist(lapply(model.list, class))=="houwie")]
  }
  if(class(model.list) != "list" | length(model.list) < 2){
    stop("getModelAvgParams requires multiple houwie model objects to be input as a list.", call. = FALSE)
  }
  rate_cats <- simplify2array(lapply(model.list, "[[", "rate.cat"))
  n_states <- simplify2array(lapply(model.list, function(x) dim(x$index.disc)[1]))
  n_obs <- unique(n_states/rate_cats)
  
  # name the models
  if(is.null(names(model.list))){
    mod_names <- paste0("M", 1:length(model.list))
    names(model.list) <- mod_names
  }else{
    mod_names <- names(model.list)
  }
  
  # pull the aic weights
  mods_table <- getModelTable(model.list)
  if(diff(range(mods_table$AIC)) > 1e5){
    if(!force){
      max_aic <- max(mods_table$AIC)
      model.list <- model.list[abs(mods_table$AIC - max_aic)  < 1e5]
      mods_table <- getModelTable(model.list)
      mod_names <- names(model.list)
    }else{
      warning("It is possible that one or more of your models failed to converge. The AIC between the best and worst models exceeds 1e10. Set force=FALSE to automatically remove potentially failed runs.")
    }
  }
  AICwts <- mods_table$AICwt
  tip_values_by_model <- lapply(model.list, get_tip_values)
  for(i in 1:length(tip_values_by_model)){
    tip_values_by_model[[i]] <- tip_values_by_model[[i]] * AICwts[i]
  }
  weighted_tip_values <- Reduce("+", tip_values_by_model)
  observed_tip_states <- model.list[[1]]$hOUwie.dat$PossibleTraits[as.numeric(model.list[[1]]$hOUwie.dat$data.cor[,2])]
  names(observed_tip_states) <- model.list[[1]]$hOUwie.dat$data.cor[,1]
  weighted_tip_values <- weighted_tip_values[match(names(observed_tip_states), rownames(weighted_tip_values)),]
  weighted_tip_values$tip_state <- observed_tip_states
  return(weighted_tip_values)
}

# different OU models have different parameter structures. This will evaluate the appropriate one.
getOUParamStructure <- function(model, nObsState, rate.cat=1, null.model=FALSE){
  if(null.model & rate.cat==1){
    cat("\nNull model was set to be true, but rate category was set to be 1. Rate category number has been set to 2.\n")
    rate.cat <- 2
  }
  if(rate.cat > 1){
    StateNames <- paste("(", rep(1:nObsState, rate.cat), rep(LETTERS[1:rate.cat], each = nObsState), ")", sep = "")
  }else{
    StateNames <- paste("(", rep(1:nObsState, rate.cat), ")", sep = "")
  }
  if(null.model){
    nState <- rate.cat
  }else{
    nState <- nObsState * rate.cat
  }
  index.mat <- matrix(NA, 3, nState, dimnames = list(c("alpha", "sigma2", "theta")))
  if(model == "BM1"){
    index.mat[2,] <- 1
    index.mat[3,] <- 2
  }
  if(model == "BMV"){
    index.mat[2,] <- 1:nState
    index.mat[3,] <- nState+1
  }
  if(model == "OU1"){
    index.mat[1,] <- 1
    index.mat[2,] <- 2
    index.mat[3,] <- 3
  }
  if(model == "OUA"){
    index.mat[1,] <- 1:nState
    index.mat[2,] <- nState+1
    index.mat[3,] <- nState+2
  }
  if(model == "OUV"){
    index.mat[1,] <- 1
    index.mat[2,] <- 2:(nState+1)
    index.mat[3,] <- nState+2
  }
  if(model == "OUM"){
    index.mat[1,] <- 1
    index.mat[2,] <- 2
    index.mat[3,] <- 3:(nState+2)
  }
  if(model == "OUVA") {
    index.mat[1,] <- 1:nState
    index.mat[2,] <- (nState+1):(nState+nState)
    index.mat[3,] <- (nState+nState)+1
  }
  if(model == "OUMA") {
    index.mat[1,] <- 1:nState
    index.mat[2,] <- (nState+1)
    index.mat[3,] <- (nState+2):(nState+nState+1)
  }
  if(model == "OUMV") {
    index.mat[1,] <- 1
    index.mat[2,] <- 2:(nState+1)
    index.mat[3,] <- (nState+2):(nState+nState+1)
  }
  if(model == "OUMVA") {
    index.mat[1,] <- 1:nState
    index.mat[2,] <- (nState+1):(nState+nState)
    index.mat[3,] <- (nState+nState+1):(nState+nState+nState)
  }
  if(null.model){
    index.mat <- index.mat[,rep(1:rate.cat, each = nObsState)]
  }
  colnames(index.mat) <- StateNames 
  return(index.mat)
}


### print function
# print a houwie object
print.houwie <- function(x, ...){
  ntips <- Ntip(x$phy)
  output <- data.frame(x$loglik,x$DiscLik, x$ContLik, x$AIC,x$AICc,x$BIC, ntips, x$param.count, row.names="")
  names(output) <- c("lnLTot","lnLDisc", "lnLCont", "AIC","AICc","BIC","nTaxa","nPars")
  cat("\nFit\n")
  print(output)
  cat("\nLegend\n")
  print(x$legend)
  cat("\nRegime Rate matrix\n")
  print(x$solution.disc)
  cat("\nOU Estimates\n")
  print(x$solution.cont)
  cat("\n")
  if(!any(is.na(x$solution.cont[1,]))){
    cat("\nHalf-life (another way of reporting alpha)\n")
    print(log(2)/x$solution.cont[1,])
  }
}

silence <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

