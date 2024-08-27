##### Main internal functions ##### 
hOUwie.dev <- function(p, phy, data, rate.cat, tip.fog,
                       index.disc, index.cont, root.p,
                       edge_liks_list, nSim, all.paths=NULL, 
                       sample_tips=FALSE, sample_nodes=FALSE,
                       adaptive_sampling=FALSE, split.liks=FALSE, 
                       global_liks_mat=global_liks_mat, diagn_msg=FALSE){
  tip.paths <- all.paths[1:length(phy$tip.label)]
  p <- exp(p)
  # check if these parameters exist in the global matrix
  # set(global_liks_mat, i = as.integer(1),  j = 1:4, value=as.list(c(0, p)))
  if(!is.null(global_liks_mat)){
    liks_match_vector <- colSums(t(global_liks_mat[,-1]) - p) == 0
    llik_houwie <- as.numeric(global_liks_mat[which(liks_match_vector), 1])
    if(!split.liks){
      if(any(liks_match_vector, na.rm = TRUE)){
        # print(llik_houwie)
        # print(p)
        return(-llik_houwie)
      }
    }
  }
  
  k <- max(index.disc, na.rm = TRUE)
  p.mk <- p[1:k]
  p.ou <- p[(k+1):length(p)] 
  Rate.mat <- matrix(1, 3, dim(index.disc)[2])
  alpha.na <- is.na(index.cont[1,])
  index.cont[is.na(index.cont)] <- max(index.cont, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[index.cont]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  rate <- index.disc
  rate[is.na(rate)] <- k + 1
  Q <- matrix(0, dim(rate)[1], dim(rate)[2])
  Q[] <- c(p.mk, 0)[rate]
  diag(Q) <- -rowSums(Q)
  edge_liks_list_init <- edge_liks_list
  # altering the conditional probabilities based on jointly sampled decendent species
  if(sample_nodes){
    edge_liks_list <- try(getCherryConditionals(phy, data, Rate.mat, Q, edge_liks_list_init, tip.paths))
    if(inherits(edge_liks_list, what="try-error")){
    #if(class(edge_liks_list) == "try-error"){
      return(1e10)
    }
  }
  # a way to alter the conditional values of the tips based on current parameeter values (only meaningful if there are multiple rate cats)
  if(rate.cat > 1 & sample_tips){
    normal.params <- rbind(theta, sigma.sq)
    sample.tip.probs <- apply(normal.params, 2, function(x) dnorm(data[,3], x[1], sqrt(x[2])))
    for(i in 1:length(phy$tip.label)){
      if(all(sample.tip.probs[i,] == 0)){
        sample.tip.probs[i,] <- rep(1, length(sample.tip.probs[i,]))
      }
      sample_tip_i <- edge_liks_list[[i]][1,] * sample.tip.probs[i,]
      edge_liks_list[[i]][1,] <- sample_tip_i/sum(sample_tip_i)
    }
  }
  # get the condtional probabilities based on the discrete values
  for(recon_index in 1:length(edge_liks_list)){
    edge_liks_list[[recon_index]] <- edge_liks_list[[recon_index]] * edge_liks_list_init[[recon_index]]
  }
  conditional_probs <- getConditionalInternodeLik(phy, Q, edge_liks_list)
  root_liks <- getRootLiks(conditional_probs, Q, root.p)
  if(is.null(root_liks)){
    return(1e10)
  }
  # initial sample
  # sample mappings based on the conditional probabilites (also calculating some time saving probabilities from transitions to and from particular states)
  internode_maps_and_discrete_probs <- getInternodeMap(phy, Q, conditional_probs$edge_liks_list, conditional_probs$root_state, root_liks, nSim, check_vector = NA, max.attempts=nSim*2)
  internode_maps <- internode_maps_and_discrete_probs$maps
  internode_samples <- internode_maps_and_discrete_probs$state_samples
  check_vector <- unlist(lapply(internode_samples, function(x) paste0(unlist(x), collapse="")))
  # if additional samples are needed (didn't reach nSim), they are taken ~randomly
  if(length(internode_samples) < nSim){
    additional_sims <- nSim - length(internode_samples)
    random_internode_maps_and_discrete_probs <- getInternodeMap(phy, Q * 100, edge_liks_list_init, conditional_probs$root_state, root_liks, additional_sims, check_vector = check_vector, max.attempts=nSim*10)
    internode_maps <- c(internode_maps, random_internode_maps_and_discrete_probs$maps)
    internode_samples <- c(internode_samples, random_internode_maps_and_discrete_probs$state_samples)
    check_vector <- unlist(lapply(internode_samples, function(x) paste0(unlist(x), collapse="")))
  }
  # calculte the discrete probabilities based on the given Q matrix (Pij already calculated)
  discrete_probs <- lapply(internode_samples, function(x) getStateSampleProb(state_sample = x, Pij = internode_maps_and_discrete_probs$Pij, root_liks = root_liks, root_edges = internode_maps_and_discrete_probs$root_edges))
  llik_discrete <- unlist(discrete_probs)
  failed_maps <- discrete_probs == -Inf
  llik_discrete <- llik_discrete[!failed_maps]
  if(length(llik_discrete) == 0){
    return(1e10)
  }
  # generate maps
  simmaps <- getMapFromSubstHistory(internode_maps[!failed_maps], phy)
  # simmaps <- lapply(simmaps, correct_map_edges)
  # if there is no character dependence the map has no influence on continuous likleihood
  character_dependence_check <- all(apply(index.cont, 1, function(x) length(unique(x)) == 1))
  if(character_dependence_check){
    llik_continuous <- OUwie.basic(simmaps[[1]], data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, tip.fog=tip.fog)
    llik_continuous <- rep(llik_continuous, length(simmaps))
  }else{
    llik_continuous <- unlist(lapply(simmaps, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, tip.fog=tip.fog)))
  }
  # combine probabilities being careful to avoid underflow
  llik_houwies <- llik_discrete + llik_continuous
  llik_houwie <- max(llik_houwies) + log(sum(exp(llik_houwies - max(llik_houwies))))
  if(!is.null(global_liks_mat)){
    if(diff(abs(c(llik_houwie,as.numeric(global_liks_mat[1,1])))) > 1e10){
      # houwie sometimes gets stuck optimizing models which have likely failed. so a quick LRT to check is implemented here
      return(1e10)
    }
    if(diff(abs(c(max(llik_discrete), max(llik_continuous))))  > 1e10){
      # another likely optimization error if the difference between the best map for discrete and continuous are over 100 orders of magnitude
      return(1e10)
    }
  }
  # after calculating the likelihoods of an intial set of maps, we sample potentially good maps
  if(adaptive_sampling & !character_dependence_check){
    adaptive_criteria <- FALSE
    adaptive_count <- 0
    best_mapping <- simmaps[[which.max(llik_houwies)]]
    while(!adaptive_criteria){
      adaptive_count <- adaptive_count + 1
      if(adaptive_count > 5){
        adaptive_criteria <- TRUE
      }
      # while we wait to meet some criteria
      # generate a set of expectations based on the best mapping
      current_ou_expectations <- getOUExpectations(best_mapping, Rate.mat, all.paths)
      # generate a new conditional probability based on the new expected values
      edge_liks_list <- try(getAdaptiveConditionals(phy, data, Rate.mat, Q, edge_liks_list_init, tip.paths, current_ou_expectations))
      if(inherits(edge_liks_list, what="try-error")){
      #if(class(edge_liks_list) == "try-error"){
        next
      }
      for(recon_index in 1:length(edge_liks_list)){
        edge_liks_list[[recon_index]] <- edge_liks_list[[recon_index]] * edge_liks_list_init[[recon_index]]
      }
      conditional_probs <- getConditionalInternodeLik(phy, Q, edge_liks_list)
      root_liks <- getRootLiks(conditional_probs, Q, root.p)
      if(is.null(root_liks)){
        next
      }
      # generate a new set of unique mappings based on the new conditional probabilities
      internode_maps_and_discrete_probs <- getInternodeMap(phy, Q, conditional_probs$edge_liks_list, conditional_probs$root_state, root_liks, nSim, check_vector = check_vector, max.attempts=nSim*2)
      if(length(internode_maps_and_discrete_probs$maps) == 0){
        # if no new maps were generated: generate the maps randomly
        internode_maps_and_discrete_probs <- getInternodeMap(phy, Q*100, edge_liks_list_init, conditional_probs$root_state, root_liks, nSim, check_vector = check_vector, max.attempts=nSim*10)
      }
      check_vector <- c(check_vector, unlist(lapply(internode_maps_and_discrete_probs$state_samples, function(x) paste0(unlist(x), collapse=""))))
      if(length(internode_maps_and_discrete_probs$maps) > 0){
        new_simmaps <- getMapFromSubstHistory(internode_maps_and_discrete_probs$maps, phy)
        if(length(new_simmaps) == 1){
          simmaps <- c(simmaps, new_simmaps[[1]])
        }else{
          simmaps <- c(simmaps, new_simmaps)
        }
        # evaluate the new mappings' joint likelhood
        discrete_probs <- lapply(internode_maps_and_discrete_probs$state_samples, function(x) getStateSampleProb(state_sample = x, Pij = internode_maps_and_discrete_probs$Pij, root_liks = root_liks, root_edges = internode_maps_and_discrete_probs$root_edges))
        continuous_probs <- lapply(new_simmaps, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, tip.fog=tip.fog))
        new_liks <- unlist(continuous_probs) + unlist(discrete_probs)
        adaptive_criteria <- max(llik_houwies) > max(new_liks)
        llik_continuous <- c(llik_continuous, unlist(continuous_probs))
        llik_discrete <- c(llik_discrete, unlist(discrete_probs))
        llik_houwies <- llik_discrete + llik_continuous
      }else{
        adaptive_criteria <- TRUE
      }
      # return to the initial step if some criteria has not been met, else done.
    }
  }
  # find the best nSim mappings after adaptive sampling
  sorted_likelihoods <- sort(llik_houwies, decreasing = TRUE, index.return = TRUE)
  unsorted_lliks_df <- data.frame(llik_discrete=llik_discrete, llik_continuous=llik_continuous)
  llik_houwies <- llik_houwies[c(na.omit(sorted_likelihoods$ix[1:nSim]))]
  llik_discrete <- llik_discrete[c(na.omit(sorted_likelihoods$ix[1:nSim]))]
  llik_continuous <- llik_continuous[c(na.omit(sorted_likelihoods$ix[1:nSim]))]
  simmaps <- simmaps[c(na.omit(sorted_likelihoods$ix[1:nSim]))]
  # get the summed probabilities using tricks to prevent underflow
  llik_houwie <- max(llik_houwies) + log(sum(exp(llik_houwies - max(llik_houwies))))
  llik_discrete_summed <- max(llik_discrete) + log(sum(exp(llik_discrete - max(llik_discrete))))
  llik_continuous_summed <- max(llik_continuous) + log(sum(exp(llik_continuous - max(llik_continuous))))
  if(!is.null(global_liks_mat)){
    if(diff(abs(c(llik_houwie,as.numeric(global_liks_mat[1,1])))) > 1e10){
      # houwie sometimes gets stuck optimizing models which have likely failed. so a quick LRT to check is implemented here
      return(1e10)
    }
  }
  if(split.liks){
    # expected_vals <- lapply(simmaps, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, tip.fog=tip.fog,return.expected.vals=TRUE))
    # expected_vals <- colSums(do.call(rbind, expected_vals) * exp(llik_houwies - max(llik_houwies))/sum(exp(llik_houwies - max(llik_houwies))))
    if(!is.na(as.numeric(global_liks_mat[which(liks_match_vector), 1]))){
      llik_houwie <- as.numeric(global_liks_mat[which(liks_match_vector), 1])
    }
    return(list(TotalLik = llik_houwie, DiscLik = llik_discrete_summed, ContLik = llik_continuous_summed, llik_discrete=llik_discrete, llik_continuous=llik_continuous, simmaps=simmaps, unsorted_lliks_df=unsorted_lliks_df))
  }
  if(!is.null(global_liks_mat)){
    new_row <- which(global_liks_mat$X1 == 0)[1]
    set(global_liks_mat, as.integer(new_row), names(global_liks_mat), as.list(c(llik_houwie, p)))
  }
  if(diagn_msg){
    print(c(round(llik_houwie, 2), round(llik_discrete_summed, 2), round(llik_continuous_summed, 2), round(p, 2)))
  }
  # print(p)
  return(-llik_houwie)
}

hOUwie.fixed.dev <- function(p, simmaps, data, rate.cat, tip.fog,
                             index.disc, index.cont, root.p, 
                             edge_liks_list, all.paths=NULL, 
                             sample_tips=FALSE, sample_nodes=FALSE,
                             split.liks=FALSE, adaptive_sampling=FALSE,
                             global_liks_mat=global_liks_mat, diagn_msg=FALSE){
  tip.paths <- all.paths[1:length(simmaps[[1]]$tip.label)]
  p <- exp(p)
  # check if these parameters exist in the global matrix
  # set(global_liks_mat, i = as.integer(1),  j = 1:4, value=as.list(c(0, p)))
  if(!is.null(global_liks_mat)){
    liks_match_vector <- colSums(t(global_liks_mat[,-1]) - p) == 0
    llik_houwie <- as.numeric(global_liks_mat[which(liks_match_vector), 1])
    if(!split.liks){
      if(any(liks_match_vector, na.rm = TRUE)){
        # print(llik_houwie)
        # print(p)
        return(-llik_houwie)
      }
    }
  }
  
  k <- max(index.disc, na.rm = TRUE)
  p.mk <- p[1:k]
  p.ou <- p[(k+1):length(p)] 
  Rate.mat <- matrix(1, 3, dim(index.disc)[2])
  alpha.na <- is.na(index.cont[1,])
  index.cont[is.na(index.cont)] <- max(index.cont, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[index.cont]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  rate <- index.disc
  rate[is.na(rate)] <- k + 1
  Q <- matrix(0, dim(rate)[1], dim(rate)[2])
  Q[] <- c(p.mk, 0)[rate]
  diag(Q) <- -rowSums(Q)
  # calculte the discrete probabilities based on the given Q matrix (Pij already calculated)
  if(inherits(root.p[1], what="character")){
  #if(class(root.p)[1] == "character"){
    if(root.p == "yang"){
      root_liks <- c(MASS::Null(Q))
      root_liks <- root_liks/sum(root_liks)
    }
    if(root.p == "flat"){
      root_liks <- rep(1/dim(Q)[1], dim(Q)[1])
    }
  }else{
    root_liks <- root.p/sum(root.p)
  }
  
  discrete_probs <- lapply(simmaps, function(x) getMapProb(x, Q, root_liks))
  llik_discrete <- unlist(discrete_probs)
  failed_maps <- discrete_probs == -Inf
  llik_discrete <- llik_discrete[!failed_maps]
  if(length(llik_discrete) == 0){
    return(1e10)
  }
  # if there is no character dependence the map has no influence on continuous likleihood
  character_dependence_check <- all(apply(index.cont, 1, function(x) length(unique(x)) == 1))
  if(character_dependence_check){
    llik_continuous <- OUwie.basic(simmaps[[1]], data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, tip.fog=tip.fog)
    llik_continuous <- rep(llik_continuous, length(simmaps))
  }else{
    llik_continuous <- unlist(lapply(simmaps, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, tip.fog=tip.fog)))
  }
  # combine probabilities being careful to avoid underflow
  llik_houwies <- llik_discrete + llik_continuous
  llik_houwie <- max(llik_houwies) + log(sum(exp(llik_houwies - max(llik_houwies))))
  llik_discrete_summed <- max(llik_discrete) + log(sum(exp(llik_discrete - max(llik_discrete))))
  llik_continuous_summed <- max(llik_continuous) + log(sum(exp(llik_continuous - max(llik_continuous))))
  
  # after calculating the likelihoods of an intial set of maps, we sample potentially good maps
  # find the best nSim mappings after adaptive sampling
  if(split.liks){
    # expected_vals <- lapply(simmaps, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, tip.fog=tip.fog,return.expected.vals=TRUE))
    # expected_vals <- colSums(do.call(rbind, expected_vals) * exp(llik_houwies - max(llik_houwies))/sum(exp(llik_houwies - max(llik_houwies))))
    if(!is.na(as.numeric(global_liks_mat[which(liks_match_vector), 1]))){
      llik_houwie <- as.numeric(global_liks_mat[which(liks_match_vector), 1])
    }
    return(list(TotalLik = llik_houwie, DiscLik = llik_discrete_summed, ContLik = llik_continuous_summed, llik_discrete=llik_discrete, llik_continuous=llik_continuous, simmaps=simmaps))
  }
  if(!is.null(global_liks_mat)){
    new_row <- which(global_liks_mat$X1 == 0)[1]
    set(global_liks_mat, as.integer(new_row), names(global_liks_mat), as.list(c(llik_houwie, p)))
  }
  if(diagn_msg){
    print(c(round(llik_houwie, 2), round(llik_discrete_summed, 2), round(llik_continuous_summed, 2), round(p, 2)))
  }
  # print(p)
  return(-llik_houwie)
}

# internal for houwie.thorough
runSingleThorough <- function(houwie_obj, new_maps, init_pars){
  hOUwie.dat <- houwie_obj$hOUwie.dat
  root.p <- houwie_obj$root.p
  tip.fog <- houwie_obj$tip.fog
  rate.cat <- houwie_obj$rate.cat
  index.disc <- houwie_obj$index.disc
  n_p_trans <- max(index.disc, na.rm = TRUE)
  p_disc <- na.omit(c(init_pars$mod_avg_disc))[1:n_p_trans]
  index.cont <- houwie_obj$index.cont
  n_p_alpha <- length(unique(na.omit(index.cont[1,])))
  n_p_sigma <- length(unique(na.omit(index.cont[2,])))
  n_p_theta <- length(unique(na.omit(index.cont[3,])))
  p_cont <- c(init_pars$mod_avg_cont[1,seq_len(n_p_alpha)],
              init_pars$mod_avg_cont[2,seq_len(n_p_sigma)],
              init_pars$mod_avg_cont[3,seq_len(n_p_theta)])
  ip <- c(p_disc, p_cont)
  res <- hOUwie.fixed(simmaps = new_maps, data = hOUwie.dat$data.ou, rate.cat = rate.cat, discrete_model = index.disc, continuous_model = index.cont, adaptive_sampling = FALSE, make_numeric = FALSE, ip = ip)
  return(res)
}

makeMapEdgesNumeric <- function(simmap, observed_traits){
  simmap$maps <- lapply(simmap$maps, function(x) replaceName(x, observed_traits))
  colnames(simmap$mapped.edge) <- as.numeric(names(observed_traits)[match(observed_traits, colnames(simmap$mapped.edge))])
  return(simmap)
}

replaceName <- function(edge, observed_traits){
  names(edge) <- as.numeric(names(observed_traits)[match(names(edge), observed_traits)])
  return(edge)
}

getEdgeLiks <- function(phy, data, n.traits, rate.cat, time_slice){
  edge_liks_list <- vector("list", dim(phy$edge)[1])
  nTip <- length(phy$tip.label)
  for(edge_i in 1:dim(phy$edge)[1]){
    # +2 because we slice the middle of the branch and need 2 terminal nodes (ancestor and descendent)
    n_slice <- (phy$edge.length[edge_i] %/% time_slice) + 2
    edge_liks_list[[edge_i]] <- matrix(1, n_slice, n.traits * rate.cat)
    if(phy$edge[edge_i,2] <= nTip){
      tmp <- numeric(n.traits)
      species_i <- phy$tip.label[phy$edge[edge_i,2]]
      state_i <- data[data[,1] == species_i, 2]
      state_i_index<- as.numeric(unlist(strsplit(as.character(state_i), "&")))
      tmp[state_i_index] <- 1
      edge_liks_list[[edge_i]][1,] <- rep(tmp, rate.cat)
    }
  }
  return(edge_liks_list)
}

getConditionalInternodeLik <- function(phy, Q, edge_liks_list){
  nTip <- length(phy$tip.label)
  external_index <- which(phy$edge[,2] <= nTip)
  # external edges
  for(edge_i in external_index){
    # move rootward along all tips to their root
    n_edges <- dim(edge_liks_list[[edge_i]])[1] - 1
    time_edge <- phy$edge.length[edge_i]
    p_mat_i <- expm(Q * (time_edge/n_edges), method=c("Ward77"))
    for(inter_edge_i in 2:(n_edges+1)){
      dec_states <- edge_liks_list[[edge_i]][inter_edge_i-1,]
      v <- edge_liks_list[[edge_i]][inter_edge_i,] * c(p_mat_i %*% dec_states )
      edge_liks_list[[edge_i]][inter_edge_i,] <- v/sum(v)
    }
  }
  # internal edges
  anc <- unique(phy$edge[,1])
  # remove the root
  root <- anc[length(anc)]
  anc <- anc[-length(anc)]
  for(anc_i in anc){
    # for the start of an internal node, combine the decs
    edge_i <- which(phy$edge[,2] == anc_i)
    dec_combo_index_i <- which(phy$edge[,1] == anc_i)
    v <- 1
    for(j in dec_combo_index_i){
      liks_j <- edge_liks_list[[j]]
      v <- v * liks_j[dim(liks_j)[1],]
    }
    v <- edge_liks_list[[edge_i]][1,] * v
    edge_liks_list[[edge_i]][1,] <- v/sum(v)
    n_edges <- dim(edge_liks_list[[edge_i]])[1] - 1
    time_edge <- phy$edge.length[edge_i]
    p_mat_i <- expm(Q * (time_edge/n_edges), method=c("Ward77"))
    for(inter_edge_i in 2:(n_edges+1)){
      dec_states <- edge_liks_list[[edge_i]][inter_edge_i-1,]
      v <- edge_liks_list[[edge_i]][inter_edge_i,] * c(p_mat_i %*% dec_states )
      edge_liks_list[[edge_i]][inter_edge_i,] <- v/sum(v)
    }
  }
  # do the root
  dec_combo_index_i <- which(phy$edge[,1] == root)
  v <- 1
  for(j in dec_combo_index_i){
    liks_j <- edge_liks_list[[j]]
    v <- v * edge_liks_list[[j]][dim(liks_j)[1],]
  }
  root_state <- v/sum(v)
  return(list(root_state = root_state,
              edge_liks_list = edge_liks_list))
}

getInternodeMap <- function(phy, Q, edge_liks_list, root_state, root_liks, nSim, check_vector=NULL, max.attempts){
  # set-up
  current.attempts <- 0
  nStates <- dim(Q)[1]
  nTip <- length(phy$tip.label)
  # a potential speedup is to calculate all Pij (bollback eq.3) for all branches first
  Pij <- array(0, c(dim(Q)[1], dim(Q)[2], length(phy$edge.length)))
  # reduced edge.lengths since we are including internodes
  number_of_nodes_per_edge <- unlist(lapply(edge_liks_list, function(x) dim(x)[1]))
  number_of_edges_per_edge <- number_of_nodes_per_edge - 1
  reduced_edge_length <- phy$edge.length/number_of_edges_per_edge
  for(i in 1:length(phy$edge.length)){
    Pij[,,i] <- expm(Q * reduced_edge_length[i])
  }
  # the probability of a descendent being in state j given starting in the row of the Pj matrix
  Pj <- vector("list", length(phy$edge.length))
  for(i in 1:length(Pj)){
    Pj[[i]] <- array(0, c(dim(Q)[1], dim(Q)[2], number_of_nodes_per_edge[i]))
    for(j in 1:number_of_nodes_per_edge[i]){
      Pj[[i]][,,j] <- sweep(Pij[,,i], MARGIN = 2, edge_liks_list[[i]][j,], '*') 
    }
  }
  # simulate nSim substitution histories
  rev.pruning.order <- rev(reorder.phylo(phy, "pruningwise", index.only = TRUE))
  sub_histories <- vector("list", nSim)
  root_edges <- which(phy$edge[,1] == nTip + 1)
  edge_index <- phy$edge
  Map_i <- mapply(function(x, y) rep(x, y), x=reduced_edge_length/2, y=number_of_edges_per_edge*2, SIMPLIFY = FALSE)
  if(!is.null(check_vector)){
    state_samples <- vector("list", nSim)
    sim_counter <- 0
    while(!(sim_counter >= nSim | current.attempts >= max.attempts)){
      state_sample <- try(getInternodeStateSample(Pj, root_state, root_edges, rev.pruning.order, edge_index, nStates, number_of_nodes_per_edge), silent = TRUE)
      if(class(state_sample) == "try-error"){
        current.attempts <- current.attempts + 1
      }else{
        current_mapping_id <- paste0(unlist(state_sample), collapse="")
        if(!current_mapping_id %in% check_vector){
          sim_counter <- sim_counter + 1
          state_samples[[sim_counter]] <- state_sample
          check_vector <- c(check_vector, current_mapping_id)
        }else{
          current.attempts <- current.attempts + 1
        }
      }
    }
  }else{
    state_samples <- lapply(1:nSim, function(x) try(getInternodeStateSample(Pj, root_state, root_edges, rev.pruning.order, edge_index, nStates, number_of_nodes_per_edge), silent = TRUE))
    state_samples <- state_samples[unlist(lapply(state_samples, class)) != "try-error"]
  }
  mapping_ids <- unlist(lapply(state_samples, function(x) paste0(unlist(x), collapse="")))
  state_samples <- state_samples[!duplicated(mapping_ids, nmax = 1)]
  mapping_ids <- unlist(lapply(state_samples, function(x) paste0(unlist(x), collapse="")))
  state_samples <- state_samples[!mapping_ids == ""]
  maps <- lapply(state_samples, function(x) getMapFromStateSample(Map_i, x))
  return(list(state_samples=state_samples, maps = maps, root_edges=root_edges,
              Pij = Pij, Pj = Pj))
}

getMapFromStateSample <- function(map, state_sample){
  for(edge_i in 1:length(map)){
    state_transitions <- rep(state_sample[[edge_i]][-c(1, length(state_sample[[edge_i]]))], each = 2)
    state_samples_i <- c(state_sample[[edge_i]][1],state_transitions,state_sample[[edge_i]][length(state_sample[[edge_i]])])
    names(map[[edge_i]]) <- state_samples_i
  }
  return(map)
}

getInternodeStateSample <- function(Pj, root_state, root_edge, rev.pruning.order, edge_index, nStates, number_of_nodes_per_edge){
  # each map will have edges split into equal time portions
  state_samples <- lapply(number_of_nodes_per_edge, function(x) numeric(x))
  root_sample <- sample(1:nStates, 1, prob = root_state)
  for(i in root_edge){
    state_samples[[i]][1] <- root_sample
  }
  # sample the nodes along a branch the last dec node goes into the next map
  for(edge_i in rev.pruning.order){
    from <- state_samples[[edge_i]][1]
    count <- 2
    n_inter_nodes <- length(state_samples[[edge_i]])
    for(inter_edge_i in (n_inter_nodes-1):1){
      to <- sample(1:nStates, 1, prob = Pj[[edge_i]][from,,inter_edge_i])
      from <- state_samples[[edge_i]][count] <- to
      count <- count + 1
    }
    ancestor_to_add <- edge_index[edge_i,2]
    anc_edge <- which(ancestor_to_add == edge_index[,1])
    for(i in anc_edge){
      state_samples[[i]][1] <- to
    }
  }
  return(state_samples)
}

# get path probability internal
getPathStateProb <- function(path_states, p_mat){
  P <- vector("numeric", length(path_states)-1)
  for(i in 1:(length(path_states)-1)){
    P[i] <- p_mat[path_states[1],path_states[2]]
    path_states <- path_states[-1]
  }
  return(sum(log(P)))
}

getStateSampleProb <- function(state_sample, Pij, root_liks, root_edges){
  path_probs <- numeric(length(state_sample))
  root_sample <- state_sample[[root_edges[1]]][1] # the root sample
  for(i in 1:length(state_sample)){
    path_probs[i] <- getPathStateProb(state_sample[[i]], Pij[,,i])
  }
  llik <- sum(path_probs) + log(root_liks[root_sample])
  return(llik)
}

# probability of a particular stochastic map
getPathMapProb <- function(map_edge, Q){
  path_states <- as.numeric(c(names(map_edge)[1], names(map_edge)[length(map_edge)]))
  P <- expm(Q * sum(map_edge))[path_states[1],path_states[2]]
  return(log(P))
}

getMapProb <- function(simmap, Q, root_prior){
  path_probs <- lapply(simmap$maps, function(x) getPathMapProb(x, Q))
  root_state <- as.numeric(names(simmap$maps[[which.min(simmap$edge[,1])]][1]))
  p_vec <- unlist(path_probs)
  llik <- sum(c(p_vec, log(root_prior)[root_state]))
  # llik <- sum(pathway_liks)
  return(llik)
}



# take substition histories and make them simmaps
getMapFromSubstHistory <- function(maps, phy){
  mapped.edge <- lapply(maps, function(x) corHMM:::convertSubHistoryToEdge(phy, x))
  obj <- vector("list", length(maps))
  for (i in 1:length(maps)){
    tree.simmap <- phy
    tree.simmap$maps <- maps[[i]]
    tree.simmap$mapped.edge <- mapped.edge[[i]]
    attr(tree.simmap, "map.order") <- "right-to-left"
    if (!inherits(tree.simmap, "simmap")) 
      class(tree.simmap) <- c("simmap", setdiff(class(tree.simmap), 
                                                "simmap"))
    obj[[i]] <- tree.simmap
  }
  if (length(maps) > 1) {
    class(obj) <- c("multiSimmap", "multiPhylo")
  }
  return(obj)
  
}
# a basic optimization for OUwie basic
OUwie.basic.dev <- function(p, phy, data, tip.fog, index.cont, tip.paths=NULL){
  p <- exp(p)
  index.cont[is.na(index.cont)] <- max(index.cont, na.rm = TRUE) + 1
  Rate.mat <- matrix(1, 3, dim(index.cont)[2])
  Rate.mat[] <- c(p, 1e-10)[index.cont]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  llik_continuous <- OUwie.basic(phy, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, tip.fog=tip.fog)
  return(-llik_continuous)
}


# probability of the continuous parameter
OUwie.basic<-function(phy, data, simmap.tree=TRUE, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, get.root.theta=FALSE, shift.point=0.5, alpha, sigma.sq, theta, tip.fog="none", algorithm="three.point", tip.paths=NULL, return.expected.vals=FALSE){
  # organize tip states based on what the simmap suggests
  mapping <- unlist(lapply(phy$maps, function(x) names(x[length(x)])))
  nTip <- length(phy$tip.label)
  TipStates <- mapping[match(match(data[,1], phy$tip.label), phy$edge[,2])]
  data[,2] <- TipStates
  
  #Makes sure the data is in the same order as the tip labels
  if(tip.fog=="none"){
    data <- data.frame(data[,2], data[,3], row.names=data[,1])
    data <- data[phy$tip.label,]
  }
  if(tip.fog=="known"){
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
  root.edge.index <- which(phy$edge[,1] == ntips+1)
  root.state <- which(colnames(phy$mapped.edge)==names(phy$maps[[root.edge.index[2]]][1]))
  ##Begins the construction of the edges matrix -- similar to the ouch format##
  edges <- cbind(c(1:(n-1)),phy$edge,MakeAgeTable(phy, root.age=root.age))
  if(scaleHeight == TRUE){
    Tmax <- max(MakeAgeTable(phy, root.age=root.age))
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
  
  map <- phy$maps
  
  Rate.mat <- rbind(alpha, sigma.sq, theta)
  pars <- matrix(c(theta, sigma.sq, alpha), length(theta), 3, dimnames = list(1:length(sigma.sq), c("opt", "sig", "alp")))
  # if the simmap did not simulate every possible state in a given hmm
  if(dim(phy$mapped.edge)[2] != dim(Rate.mat)[2]){
    Rate.mat <- Rate.mat[,as.numeric(colnames(phy$mapped.edge))]
    pars <- pars[as.numeric(colnames(phy$mapped.edge)), ]
  }
  
  if(get.root.theta == TRUE){
    W <- weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=FALSE, shift.point=shift.point)
  }else{
    W <- weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=TRUE, shift.point=shift.point)
  }
  
  #Likelihood function for estimating model parameters
  if(get.root.theta == TRUE){
    root_index <- as.numeric(names(phy$maps[[which.min(phy$edge[,1])[1]]]))
    theta0 <- theta[root_index]
    expected.vals <- colSums(t(W) * c(theta0, pars[,1]))
    names(expected.vals) <- phy$tip.label
  }else{
    expected.vals <- colSums(t(W) * pars[,1])
    names(expected.vals) <- phy$tip.label
  }
  if(return.expected.vals){
    return(expected.vals)
  }
  transformed.tree <- transformPhy(phy, map, pars, tip.paths)
  # generate a map from node based reconstructions
  if(tip.fog=="known"){
    TIPS <- transformed.tree$tree$edge[,2] <= length(transformed.tree$tree$tip.label)
    transformed.tree$tree$edge.length[TIPS] <- transformed.tree$tree$edge.length[TIPS] + (data[,3]^2/transformed.tree$diag/transformed.tree$diag)
  }
  comp <- NA
  try(comp <- phylolm::three.point.compute(transformed.tree$tree, x, expected.vals, transformed.tree$diag, check.precision = FALSE), silent=TRUE) # for the newest version (github) of phylolm
  if(is.na(comp[1])){
    try(comp <- phylolm::three.point.compute(transformed.tree$tree, x, expected.vals, transformed.tree$diag), silent=TRUE) # for the cran version of phylolm
  }
  if(is.na(comp[1])){
    return(-1e10)
  }else{
    nTips <- length(phy$tip.label)
    logl <- -as.numeric(nTips * log(2 * pi) + comp$logd + (comp$PP - 2 * comp$QP + comp$QQ))/2
    return(logl)
  }
}

# script for generating all the possible underlying mappings and looking at joint probablity. using this we can look at the bias produced by looking only at the discrete mappings.
fixEdgeLiksLiks <- function(edge_liks_list, combo, phy, n_tips, n_nodes, n_internodes, nStates, rate.cat){
  # fix the externals
  for(j in 1:n_tips){
    tip_index <- n_nodes + n_internodes + which(phy$edge[phy$edge[,2] <= n_tips,2] == j)
    tip_state <- combo[tip_index]
    dec_edges_to_fix <- which(phy$edge[,2] == j)
    fix_vector <- numeric(nStates * rate.cat)
    fix_vector[tip_state] <- 1
    for(k in dec_edges_to_fix){
      edge_liks_list[[k]][1,] <- fix_vector
    }
  }
  # fix the internals
  for(j in 1:n_nodes){
    node_index <- unique(phy$edge[,1])[j]
    node_state <- combo[j]
    anc_edges_to_fix <- which(phy$edge[,1] == node_index)
    dec_edges_to_fix <- which(phy$edge[,2] == node_index)
    fix_vector <- numeric(nStates * rate.cat)
    fix_vector[node_state] <- 1
    for(k in dec_edges_to_fix){
      edge_liks_list[[k]][1,] <- fix_vector
    }
    for(k in anc_edges_to_fix){
      last_row <- dim(edge_liks_list[[k]])[1]
      edge_liks_list[[k]][last_row,] <- fix_vector
    }
  }
  # fix the inter nodes
  if(n_internodes > 0){
    for(j in 1:n_internodes){
      internode_index_list <- which(unlist(lapply(edge_liks_list, function(x) dim(x)[1] - 2)) >= 1)[j]
      internode_state <- combo[j + n_nodes]
      internode_index_edge <- which(apply(edge_liks_list[[internode_index_list]], 1, function(x) sum(x) > 1))[1]
      fix_vector <- numeric(nStates * rate.cat)
      fix_vector[internode_state] <- 1
      edge_liks_list[[internode_index_list]][internode_index_edge,] <- fix_vector
    }
  }
  return(edge_liks_list)
}

getAllJointProbs<- function(phy, data, rate.cat, time_slice, Q, alpha, sigma.sq, theta, quiet=TRUE){
  # prerequisites
  hOUwie.dat <- organizeHOUwieDat(data, "none", TRUE)
  nStates <- as.numeric(max(hOUwie.dat$data.cor[,2]))
  tip.paths <- lapply(1:length(phy$tip.label), function(x) getPathToRoot(phy, x))
  # generate the edge_liks_list
  edge_liks_list <- getEdgeLiks(phy, hOUwie.dat$data.cor, nStates, rate.cat, time_slice)
  # determine all the possible ways to fix the node states
  # how many nodes and internodes are we fixing?
  n_tips <- length(phy$tip.label)
  n_nodes <- n_tips - 1
  n_internodes <- sum(unlist(lapply(edge_liks_list, function(x) dim(x)[1] - 2)))
  # what are the possible internal combinations?
  internal_possibilities <- rep(list(1:(nStates*rate.cat)), n_nodes + n_internodes)
  external_possibilities <- lapply(edge_liks_list[phy$edge[,2] <= n_tips], function(x) which(x[1,] == 1))
  all_combinations <- expand.grid(c(internal_possibilities, external_possibilities))
  # the joint probability table
  joint_probability_table <- matrix(NA, dim(all_combinations)[1], 3, dimnames = list(1:dim(all_combinations)[1], c("disc", "cont", "total")))
  simmap_list <- vector("list", dim(all_combinations)[1])
  if(!quiet){
    cat("Begining to calcualte all possible map combinations...\n")
  }
  # for each possibility, generate an edge_liks_list to match
  for(i in 1:dim(all_combinations)[1]){
    if(!quiet){
      cat("\r", i, "of", dim(all_combinations)[1], "complete.         ")
    }
    combo_i <- as.numeric(all_combinations[i,])
    root_state <- numeric(nStates * rate.cat)
    root_state[combo_i[which(unique(phy$edge[,1]) == n_tips + 1)[1]]] <- 1
    edge_liks_list_i <- fixEdgeLiksLiks(edge_liks_list, combo_i, phy, n_tips, n_nodes, n_internodes, nStates, rate.cat)
    root_liks <- rep(1, nStates * rate.cat)/(nStates * rate.cat)
    # calculate the discrete probability of the edge_liks_list
    tmp <- getInternodeMap(phy, Q, edge_liks_list_i, root_state, root_liks, 1)
    internode_maps <- tmp$maps
    internode_samples <- tmp$state_samples
    llik_discrete <- unlist(lapply(internode_samples, function(x) getStateSampleProb(state_sample = x, Pij = tmp$Pij, root_liks = root_liks, root_edges = tmp$root_edges)))
    # generate a stochstic map
    simmaps <- getMapFromSubstHistory(internode_maps, phy)
    llik_continuous <- unlist(lapply(simmaps, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, tip.fog="none")))
    simmap_list[[i]] <- simmaps[[1]]
    joint_probability_table[i,] <- c(llik_discrete, llik_continuous, llik_discrete + llik_continuous)
  }
  if(!quiet){
    cat("\n")
  }
  simmap_list <- lapply(simmap_list, correct_map_edges)
  class(simmap_list) <- c("multiSimmap", "multiPhylo")
  return(list(joint_probability_table=joint_probability_table, simmap_list=simmap_list, all_combinations=all_combinations))
}

# branchwise joint calculations for conditional probabilities
getOUProbBranch <- function(tip_value, states, pars, bl, init_value, init_var){
  # states are read from tip to root states[1] = tip value, states[2] = node value
  times <- c(0, bl/2, bl)
  alphas <- c(pars[1,states])
  sigma2s <- c(pars[2,states])
  thetas <- c(pars[3,states])
  if(is.null(init_value)){
    values <- c(thetas[1], thetas)
  }else{
    values <- c(init_value, thetas)
  }
  # finding the expected value
  tip_weight <- exp(-((alphas[1] * (times[2] - times[1])) + (alphas[2] * (times[3] - times[2]))))
  theta_1_weight <- tip_weight * (exp(alphas[1] * times[2]) - exp(alphas[1] * times[1]))
  theta_2_weight <- tip_weight * (exp(alphas[2] * times[3]) - exp(alphas[2] * times[2]))
  weights <- c(tip_weight, theta_1_weight, theta_2_weight)/sum(c(tip_weight, theta_1_weight, theta_2_weight))
  expected_value <- sum(values * weights)
  # finding the variance 
  var_weight <- exp(-((2 * alphas[1] * (times[2] - times[1])) + (2 * alphas[2] * (times[2] - times[1]))))
  var_1 <- sigma2s[1]/(2*alphas[1]) * (exp(2 * alphas[1] * times[2]) - exp(2 * alphas[1] * times[1]))
  var_2 <- sigma2s[2]/(2*alphas[2]) * (exp(2 * alphas[2] * times[3]) - exp(2 * alphas[2] * times[2]))
  variance <- (sum(c(var_1, var_2)) * var_weight) + init_var
  loglik <- dnorm(tip_value, expected_value, sqrt(variance), log = TRUE)
  return(loglik)
}

getJointProbBranch <- function(states, tip_value, pars, bl, P_mat, init_value, init_var){
  coninuous_prob <- getOUProbBranch(tip_value, states, pars, bl, init_value, init_var)
  discrete_prob <- log(P_mat[states[1], states[2]])
  joint_prob <- sum(coninuous_prob, discrete_prob)
  return(joint_prob)
}

getJointBranchMatrix <- function(possible_internal, possible_external, tip_value, pars, bl, P_mat, init_value=NULL, init_var = 0){
  cond_matrix <- matrix(NA, dim(P_mat)[1], dim(P_mat)[2])
  for(i in possible_internal){
    for(j in possible_external){
      cond_matrix[i,j] <- getJointProbBranch(c(i,j), tip_value, pars, bl, P_mat, init_value, init_var)
    }
  }
  colnames(cond_matrix) <- rownames(cond_matrix) <- c(1:max(possible_internal))
  return(cond_matrix)
}

# a function for getting OU expectations per node and tip: both variance and means
getOUExpectations <- function(simmap, Rate.mat, tip.paths=NULL){
  Rate.mat[is.na(Rate.mat)] <- 1e-10
  # phy must be of class simmap
  nTip <- length(simmap$tip.label)
  if(is.null(tip.paths)){
    tip.paths <- lapply(1:(nTip*2-1), function(x) getPathToRoot(simmap, x))
  }
  # find the root state and its theta value
  root_state <- as.numeric(names(simmap$maps[[which.min(simmap$edge[,1])]])[1])
  expected_vars <- expected_values <- vector("numeric", length(tip.paths))
  for(i in 1:length(tip.paths)){
    # evaluate the map for this particular edge and calculate the tipward variance
    Map_i <- simmap$maps[tip.paths[[i]]]
    if(length(Map_i) == 0){
      next # we are at the root
    }
    # get the branch specific parameters
    branch_lengths_i <- unlist(rev(Map_i))
    alpha_values <- Rate.mat[1, as.numeric(names(branch_lengths_i))]
    sigma_values <- Rate.mat[2, as.numeric(names(branch_lengths_i))]
    theta_values <- Rate.mat[3, as.numeric(names(branch_lengths_i))]
    # get the ancestral weight value
    ancestral_weight <- exp(-sum(branch_lengths_i * alpha_values))
    # get the ancestral weight variance
    ancestral_var <- exp(-sum(2 * branch_lengths_i * alpha_values))
    # get the weight of each theta per branch
    dist_from_root_i <- c(0, cumsum(branch_lengths_i))
    var_weights <- theta_weights <- vector("numeric", length(dist_from_root_i))
    for(j in 2:length(dist_from_root_i)){
      # doing means
      time_2_j <- exp(dist_from_root_i[j] * alpha_values[j-1])
      time_1_j <- exp(dist_from_root_i[j-1] * alpha_values[j-1])
      weight_j <- ancestral_weight * (time_2_j - time_1_j)
      theta_weights[j] <- weight_j
      # doing variance
      var_2 <- exp(2 * dist_from_root_i[j] * alpha_values[j-1])
      var_1 <- exp(2 * dist_from_root_i[j-1] * alpha_values[j-1])
      var_weights[j] <- (sigma_values[j-1]/(2 * alpha_values[j-1])) * (var_2 - var_1)
    }
    # doing the means
    theta_weights[1] <- ancestral_weight # add the ancestral weight 
    theta_weights <- theta_weights/sum(theta_weights) # standardize the weights to sum to one 
    theta_values <- c(Rate.mat[3, root_state], theta_values) # add the root theta to the branch values
    expected_value <- sum(theta_values * theta_weights)
    expected_values[i] <- expected_value
    # doing the variance
    expected_vars[i] <- sum(var_weights) * ancestral_var
  }
  expected_values[expected_values == 0] <- Rate.mat[3, root_state] # add the root as an expecation
  names(expected_vars) <- names(expected_values) <- 1:(nTip*2-1) # named by node number
  return(list(expected_means=expected_values, expected_variances=expected_vars))
}

combineDesc <- function(P_mat_1, P_mat_2){
  # P_mat_1 <- node_mat[[1]]
  # P_mat_2 <- node_mat[[2]]
  normalized_combined_probs <- vector("list", dim(P_mat_1)[1])
  for(i in 1:dim(P_mat_1)[1]){
    combined_probs_i <- P_mat_1[i,] * t(P_mat_2)
    normalized_combined_probs[[i]] <- t(combined_probs_i)/colSums(combined_probs_i)
  }
  normalized_combined_probs <- do.call(rbind, normalized_combined_probs)
  return(normalized_combined_probs)  
}


# a function for generating an altnerative conditional probabilities (cherry sampling of continuous)
getCherryConditionals <- function(phy, data, Rate.mat, Q, edge_liks_list,tip.paths){
  node_ages <- branching.times(phy)
  dec_liks <- do.call(rbind, lapply(edge_liks_list, function(x) x[1,]))
  anc <- unique(phy$edge[,1])
  possible_internal <- 1:dim(Q)[1]
  for(i in 1:length(anc)){
    # whatever tip paths are chosen, they must include each of these
    # i = 3
    node_edges <- which(phy$edge[,1] == anc[i])
    bl <- node_ages[names(node_ages) == anc[i]]
    node_mat <- vector("list", length(node_edges))
    P_mat <- expm(Q * bl)
    for(j in 1:length(node_edges)){
      # j = 1
      # tips which contain j are shown here
      tips_from_anc <- which(unlist(lapply(tip.paths, function(x) any(node_edges[j] %in% x))))
      # tip_sampled <- sample(c(tips_from_anc, tips_from_anc), 1)
      tip_sampled <- tips_from_anc
      tip_index <- match(tip_sampled, phy$edge[,2]) # relative to the edge matrix
      possible_external <- matrix(dec_liks[tip_index,], ncol = length(possible_internal))
      # possible_external <- which(dec_liks[tip_index,] == 1)
      # tip_value <- mean(data[data$sp %in% phy$tip.label[tip_sampled],3])
      tip_values <- c(data[data$sp %in% phy$tip.label[tip_sampled],3])
      # branch_matrix <- getJointBranchMatrix(possible_internal, possible_external, tip_value, Rate.mat, bl, P_mat)
      branch_matrices <- vector("list", length(tip_values))
      for(k in 1:length(tip_values)){
        branch_matrices[[k]] <- getJointBranchMatrix(possible_internal, which(possible_external[k,] == 1), tip_values[k], Rate.mat, bl, P_mat)
      }
      # branch_matrix <- getJointBranchMatrix(possible_internal, possible_internal, tip_value, Rate.mat, bl, P_mat)
      # the likelihood that the rootward state led to the known tip ward state
      node_state_liks_list <- lapply(branch_matrices, 
                                     function(x) apply(x, 1, sum_lliks))
      node_state_liks_list <- lapply(node_state_liks_list, 
                                     function(x) exp(x - max(x))/sum(exp(x - max(x))))
      node_state_liks <- do.call(rbind, node_state_liks_list)
      node_mat[[j]] <- node_state_liks
    }
    node_mat <- lapply(node_mat, colMeans)
    new_node_mat <- do.call(rbind, node_mat)
    # new_node_mat <- node_mat[[1]]
    # for(j in 2:length(node_mat)){
    #   new_node_mat <- combineDesc(new_node_mat, node_mat[[j]])
    # }
    # node_mat contains two (or more lists) of tip samples, what needs to happen now is the combination of the two edges as parent daughters and then the combination of those pairs
    # node_state_liks <- Reduce("*", node_state_liks_list)
    # node_state_liks <- apply(branch_matrix, 1, sum_lliks)
    # node_state_probs <- exp(node_state_liks - max(node_state_liks))/sum(exp(node_state_liks - max(node_state_liks)))
    # once that node has finished calculating it's conditional probabilitity of each state, we combine the two dec tips
    # node_mat[node_mat == 0] <- 1e-10
    # node_cond_prob <- apply(node_mat, 2, prod)
    node_cond_prob <- colMeans(new_node_mat)
    # node_state_probs[node_state_probs == 0] <- 1e-10
    # node_cond_prob <- node_state_probs
    # this gets place in the edge matrix
    anc_index <- which(phy$edge[,1] %in% anc[i])
    dec_index <- which(phy$edge[,2] %in% anc[i])
    for(k in anc_index){
      # k = 5
      edge_liks_list[[k]][dim(edge_liks_list[[k]])[1],] <- node_cond_prob/sum(node_cond_prob)
    }
    if(length(dec_index) > 0){
      for(k in dec_index){
        edge_liks_list[[k]][1,] <- node_cond_prob/sum(node_cond_prob)
      }
    }
  }
  return(edge_liks_list)
}

# a function for generating an altnerative conditional probabilities (adaptive sampling of continuous)
getAdaptiveConditionals <- function(phy, data, Rate.mat,  Q, edge_liks_list, tip.paths, ou_expectations){
  node_ages <- branching.times(phy)
  dec_liks <- do.call(rbind, lapply(edge_liks_list, function(x) x[1,]))
  anc <- unique(phy$edge[,1])
  possible_internal <- 1:dim(Q)[1]
  for(i in 1:length(anc)){
    # whatever tip paths are chosen, they must include each of these
    # i = 3
    anc_value <- ou_expectations$expected_means[names(ou_expectations$expected_means) == anc[i]]
    anc_var <- ou_expectations$expected_variances[names(ou_expectations$expected_variances) == anc[i]]
    node_edges <- which(phy$edge[,1] == anc[i])
    bl <- node_ages[names(node_ages) == anc[i]]
    node_mat <- vector("list", length(node_edges))
    P_mat <- expm(Q * bl)
    for(j in 1:length(node_edges)){
      # j = 1
      # tips which contain j are shown here
      tips_from_anc <- which(unlist(lapply(tip.paths, function(x) node_edges[j] %in% x)))
      # tips_from_anc <- which(unlist(lapply(tip.paths, function(x) any(node_edges %in% x))))
      # tip_sampled <- sample(c(tips_from_anc, tips_from_anc), 1)
      tip_sampled <- tips_from_anc
      tip_index <- match(tip_sampled, phy$edge[,2]) # relative to the edge matrix
      possible_external <- matrix(dec_liks[tip_index,], ncol = length(possible_internal))
      # tip_value <- mean(data[data$sp %in% phy$tip.label[tip_sampled],3])
      tip_values <- c(data[data$sp %in% phy$tip.label[tip_sampled],3])
      # branch_matrix <- getJointBranchMatrix(possible_internal, possible_external, tip_value, Rate.mat, bl, P_mat)
      # branch_matrices <- lapply(tip_values, function(x) getJointBranchMatrix(possible_internal, possible_external, x, Rate.mat, bl, P_mat))
      branch_matrices <- vector("list", length(tip_values))
      for(k in 1:length(tip_values)){
        branch_matrices[[k]] <- getJointBranchMatrix(possible_internal, which(possible_external[k,] == 1), tip_values[k], Rate.mat, bl, P_mat, init_value = anc_value, init_var = anc_var)
      }
      # branch_matrix <- getJointBranchMatrix(possible_internal, possible_internal, tip_value, Rate.mat, bl, P_mat, init_value = anc_value, init_var = anc_var)
      # the likelihood that the rootward state led to the known tip ward state
      node_state_liks_list <- lapply(branch_matrices, 
                                     function(x) apply(x, 1, sum_lliks))
      node_state_liks_list <- lapply(node_state_liks_list, 
                                     function(x) exp(x - max(x))/sum(exp(x - max(x))))
      node_state_liks <- do.call(rbind, node_state_liks_list)
      node_mat[[j]] <- node_state_liks
    }
    new_node_mat <- node_mat[[1]]
    for(j in 2:length(node_mat)){
      new_node_mat <- combineDesc(new_node_mat, node_mat[[j]])
    }
    node_cond_prob <- colMeans(new_node_mat)
    # once that node has finished calculating it's conditional probabilitity of each state, we combine the two dec tips
    # node_mat[node_mat == 0] <- 1e-10
    # node_cond_prob <- apply(node_mat, 2, prod)
    # this gets place in the edge matrix
    anc_index <- which(phy$edge[,1] %in% anc[i])
    dec_index <- which(phy$edge[,2] %in% anc[i])
    for(k in anc_index){
      # k = 5
      edge_liks_list[[k]][dim(edge_liks_list[[k]])[1],] <- node_cond_prob/sum(node_cond_prob)
    }
    if(length(dec_index) > 0){
      for(k in dec_index){
        edge_liks_list[[k]][1,] <- node_cond_prob/sum(node_cond_prob)
      }
    }
  }
  return(edge_liks_list)
}

# get root liks from conditionals
getRootLiks <- function(conditional_probs, Q, root.p){
  if(any(is.na(conditional_probs$root_state))){
    return(NULL)
  }
  if(inherits(root.p[1], what="character")){
  #if(class(root.p)[1] == "character"){
    if(root.p == "yang"){
      root_liks <- c(MASS::Null(Q))
      root_liks <- root_liks/sum(root_liks)
    }
    if(root.p == "flat"){
      root_liks <- rep(1/dim(Q)[1], dim(Q)[1])
    }
    if(root.p == "maddfitz"){
      root_liks <- conditional_probs$root_state/sum(conditional_probs$root_state)
    }
  }else{
    root_liks <- root.p/sum(root.p)
  }
  return(root_liks)
}

##### Utility internal functions ##### 
getDiscreteModel <- function(data, model, rate.cat, dual, collapse){
  rate <- getStateMat4Dat(data, model, dual, collapse)$rate.mat
  if (rate.cat > 1) {
    StateMats <- vector("list", rate.cat)
    for (i in 1:rate.cat) {
      StateMats[[i]] <- rate
    }
    rate <- getFullMat(StateMats)
  }
  return(rate)
}


# getAllContinuousModelStructures <- function(k, type = "OU"){
#   # index.mat <- matrix(0, 3, k, dimnames = list(c("alpha", "sigma.sq", "theta"), c(1:k)))
#   # we want all unique combinations of a parameter. then we can add a single all same
#   # how many combinations are there of 1:k numbers?
#   potential_combos <- apply(partitions:::setparts(k), 2, function(x) paste(x, collapse="_"))
#   # this technically isn't all the possible alpha combinations, but for sim purposes we're fine.
#   if(type == "BM"){
#     alpha.combos <- paste(rep(0, k), collapse="_")
#     theta.combos <- paste(rep(1, k), collapse="_")
#   }
#   if(type == "OU"){
#     alpha.combos <- potential_combos
#     theta.combos <- potential_combos
#   }
#   if(type == "BMOU"){
#     if(k > 2){
#       stop("BMOU must be manually created if k > 1 atm. Sorry.")
#     }
#     # needed_numerals <- 1:((2^k)-2)
#     needed_numerals <- 1
#     alpha.combos <- apply(sapply(needed_numerals, function(x) as.numeric(intToBits(x)[1:k])), 2, function(x) paste(x, collapse="_")) # currently doesn't allow for BM mixed with OUA
#     # theta.combos <- potential_combos
#     theta.combos <- paste(rep(1, k), collapse="_")
#   }
#   sigma.sq.combos <- potential_combos
#   all_combos <- expand.grid(list(alpha.combos, sigma.sq.combos, theta.combos))
#   index_mats <- array(NA, c(3, k, dim(all_combos)[1]), dimnames = list(c("alpha", "sigma.sq", "theta"), c(1:k)))
#   for(i in 1:dim(all_combos)[1]){
#     alpha_i <- as.numeric(unlist(strsplit(as.character(all_combos[i,1]), "_")))
#     alpha_i[alpha_i == 0] <- NA
#     sigma_i <- max(c(0, alpha_i), na.rm = TRUE) + as.numeric(unlist(strsplit(as.character(all_combos[i,2]), "_")))
#     theta_i <- max(sigma_i) + as.numeric(unlist(strsplit(as.character(all_combos[i,3]), "_")))
#     index_mats[,,i] <- rbind(alpha_i, sigma_i, theta_i)
#   }
#   return(index_mats)
# }


organizeHOUwieDat <- function(data, tip.fog, collapse = TRUE){
  # return a list of corHMM data and OU data
  if(tip.fog=="known"){
    data.cor <- data[, 1:(dim(data)[2]-2)]
    data.cor <- corHMM:::corProcessData(data.cor, collapse = collapse)
    data.ou <- data.frame(sp = data[,1], 
                          reg = data.cor$corData[,2], 
                          x = data[, dim(data)[2]-1],
                          err = data[, dim(data)[2]])
  }
  if(tip.fog=="none"){
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

organizeHOUwiePars <- function(pars, index.disc, index.cont){
  k <- max(index.disc, na.rm = TRUE)
  p.mk <- pars[1:k]
  p.ou <- pars[(k+1):length(pars)] 
  # ouwie pars
  Rate.mat <- matrix(1, 3, dim(index.disc)[2])
  index.cont[is.na(index.cont)] <- max(index.cont, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, NA)[index.cont]
  rownames(Rate.mat) <- rownames(index.cont)
  rate <- index.disc
  rate[is.na(rate)] <- k + 1
  Q <- matrix(0, dim(rate)[1], dim(rate)[2])
  Q[] <- c(p.mk, NA)[rate]
  diag(Q) <- -rowSums(Q)
  if(!is.null(colnames(rate))){
    colnames(Rate.mat) <- colnames(rate)
    colnames(Q) <- rownames(Q) <- colnames(rate)
  }
  # corhmm pars 
  return(list(solution.ou = Rate.mat,
              solution.cor = Q))
}

getIP.theta <- function(x, states, index){
  ip.theta <- vector("numeric", length(unique(index)))
  for(i in 1:length(unique(index))){
    state_i <- which(unique(index)[i] == index)
    ip.theta[i] <- mean(x[states %in% state_i])
  }
  ip.theta[is.nan(ip.theta)] <- mean(x)
  return(ip.theta)
}

simCharacterHistory <- function(phy, Q, root.freqs, Q2 = NA, NoI = NA){
  #Randomly choose starting state at root using the root.values as the probability:
  root.value <- sample.int(dim(Q)[2], 1, FALSE, prob=root.freqs/sum(root.freqs))
  #Reorder the phy:
  phy <- reorder.phylo(phy, "postorder")
  ntips <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- ntips + 1 #perhaps use an accessor to get the root node id
  
  #Generate vector that contains the simulated states:
  CharacterHistory <- integer(ntips + phy$Nnode)
  CharacterHistory[ROOT] <- as.integer(root.value)
  anc <- phy$edge[, 1]
  des <- phy$edge[, 2]
  edge.length <- phy$edge.length
  diag(Q) = 0
  diag(Q) = -rowSums(Q)
  
  # setting up the alternative Q matrix at the node of interest
  if(!any(is.na(Q2))){
    diag(Q2) = 0
    diag(Q2) = -rowSums(Q2)
  }
  if(!is.na(NoI)){
    NewQDesc <- getDescendants(phy, NoI)
  }
  
  #standard simulation protocol
  if(any(is.na(Q2)) | is.na(NoI)){
    for (i in N:1) {
      p <- expm(Q * edge.length[i], method="Ward77")[CharacterHistory[anc[i]], ]
      CharacterHistory[des[i]] <- sample.int(dim(Q)[2], size = 1, FALSE, prob = p)
    }
  }
  
  # simulating a clade under a different (Q2) evolutionary model
  if(!any(is.na(Q2)) & !is.na(NoI)){
    for (i in N:1) {
      if(anc[i] %in% NewQDesc){
        p <- expm(Q2 * edge.length[i], method="Ward77")[CharacterHistory[anc[i]], ]
        CharacterHistory[des[i]] <- sample.int(dim(Q2)[2], size = 1, FALSE, prob = p)
      }else{
        p <- expm(Q * edge.length[i], method="Ward77")[CharacterHistory[anc[i]], ]
        CharacterHistory[des[i]] <- sample.int(dim(Q)[2], size = 1, FALSE, prob = p)
      }
    }
  }
  
  TipStates <-  CharacterHistory[1:ntips]
  names(TipStates) <- phy$tip.label
  NodeStates <- CharacterHistory[ROOT:(N+1)]
  names(NodeStates) <- ROOT:(N+1)
  
  res <- list(TipStates = TipStates, NodeStates = NodeStates)
  return(res)
  #return(CharacterHistory)
}
# organize the houwie output
getHouwieObj <- function(liks_houwie, pars, phy, data, hOUwie.dat, rate.cat, tip.fog, index.disc, index.cont, root.p, nSim, sample_tips, sample_nodes, adaptive_sampling, nStates, discrete_model, continuous_model, time_slice, root.station, get.root.theta,lb_discrete_model,ub_discrete_model,lb_continuous_model,ub_continuous_model,ip, opts, quiet, negative_values){
  param.count <- max(index.disc, na.rm = TRUE) + max(index.cont, na.rm = TRUE)
  nb.tip <- length(phy$tip.label)
  solution <- organizeHOUwiePars(pars=pars, index.disc=index.disc, index.cont=index.cont)
  if(rate.cat > 1){
    StateNames <- paste("(", rep(1:nStates, rate.cat), rep(LETTERS[1:rate.cat], each = nStates), ")", sep = "")
  }else{
    StateNames <- paste("(", rep(1:nStates, rate.cat), ")", sep = "")
  }
  rownames(solution$solution.cor) <- colnames(solution$solution.cor) <- StateNames
  colnames(solution$solution.ou) <- StateNames
  names(hOUwie.dat$ObservedTraits) <- 1:length(hOUwie.dat$ObservedTraits)
  if(negative_values){
    n_theta <- length(unique(index.cont[3,]))
    pars[(length(pars) - n_theta + 1):length(pars)] <- pars[(length(pars) - n_theta + 1):length(pars)]- 50
    if(tip.fog == "none"){
      data[,dim(data)[2]] <- data[,dim(data)[2]] - 50
      solution$solution.ou[3,] <- solution$solution.ou[3,] - 50
    }else{
      data[,dim(data)[2]-1] <- data[,dim(data)[2]-1] - 50
      solution$solution.ou[3,] <- solution$solution.ou[3,] - 50
    }
  }
  obj <- list(
    loglik = liks_houwie$TotalLik,
    DiscLik = liks_houwie$DiscLik,
    ContLik = liks_houwie$ContLik,
    AIC = -2*liks_houwie$TotalLik + 2*param.count,
    AICc = -2*liks_houwie$TotalLik + 2*param.count + ((2*(param.count^2) + 2*param.count)/(nb.tip - param.count - 1)),
    BIC = -2*liks_houwie$TotalLik + log(nb.tip) * param.count,
    param.count = param.count,
    solution.disc = solution$solution.cor,
    solution.cont = solution$solution.ou,
    recon = NULL,
    index.disc = index.disc,
    index.cont = index.cont,
    phy = phy,
    legend = hOUwie.dat$ObservedTraits,
    data = data, 
    hOUwie.dat = hOUwie.dat,
    rate.cat = rate.cat, 
    discrete_model=discrete_model, 
    continuous_model=continuous_model, 
    root.p=root.p, 
    time_slice=time_slice,
    root.station=root.station, 
    get.root.theta=get.root.theta, 
    tip.fog = tip.fog, 
    sample_tips = sample_tips,
    sample_nodes = sample_nodes,
    adaptive_sampling = adaptive_sampling,
    lb_discrete_model=lb_discrete_model, 
    ub_discrete_model=ub_discrete_model,
    lb_continuous_model=lb_continuous_model, 
    ub_continuous_model=ub_continuous_model,
    p=pars, 
    ip=ip,
    nSim=nSim, 
    opts=opts,
    quiet=quiet
  )
  class(obj) <- "houwie"
  return(obj)
}
# generate random starting values for multiple starts
generateMultiStarting <- function(starts, index.disc, index.cont, n_starts, lower, upper){
  n_p_trans <- max(index.disc, na.rm = TRUE)
  n_p_alpha <- length(unique(na.omit(index.cont[1,])))
  n_p_sigma <- length(unique(na.omit(index.cont[2,])))
  n_p_theta <- length(unique(na.omit(index.cont[3,])))
  multiple_starts <- vector("list", n_starts)
  multiple_starts[[1]] <- checkStartingUBLB(starts, lower, upper)
  if(n_starts > 1){
    for(i in 2:n_starts){
      multiple_starts[[i]] <- numeric(length(starts))
      for(j in 1:length(starts)){
        if(starts[j]/10 <= lower[j] | starts[j]*10 >= upper[j]){
          lb <- (lower[j]+lower[j]*0.01)*10
          ub <- (upper[j]-upper[j]*0.01)/10
        }else{
          lb <- starts[j]
          ub <- starts[j]
        }
        multiple_starts[[i]][j] <- exp(runif(1, log(lb/10), log(ub*10)))
      }
      multiple_starts[[i]] <- checkStartingUBLB(multiple_starts[[i]], lower, upper)
    }
  }
  return(multiple_starts)
}

checkStartingUBLB <- function(starts, lower, upper){
  for(i in 1:length(starts)){
    if(starts[i] < lower[i]){
      starts[i] <- lower[i] * 2
    }
    if(starts[i] > upper[i]){
      starts[i] <- upper[i] / 2
    }
  }
  return(starts)
}

# a function for correcting edge labels to be plotted when using get all joint probs
correct_map_edges <- function(simmap){
  simmap$maps <- lapply(simmap$maps, correct_edge)
  return(simmap)
}

correct_edge <- function(edge){
  edge_names <- names(edge)
  count <- 1
  edge_merge <- numeric(length(edge))
  edge_merge[1] <- count
  for(i in 2:length(edge)){
    if(edge_names[i-1] == edge_names[i]){
      edge_merge[i] <- count
    }else{
      count <- count + 1
      edge_merge[i] <- count
    }
  }
  new_edge_lengths <- vector("numeric", length(unique(edge_merge)))
  new_edge_names <- vector("character", length(unique(edge_merge)))
  for(j in unique(edge_merge)){
    new_edge_lengths[j] <- sum(edge[edge_merge %in% j])
    new_edge_names[j] <- names(edge)[edge_merge %in% j][1]
  }
  names(new_edge_lengths) <- new_edge_names
  return(new_edge_lengths)
}

# simple funciton for getting a probability from logliks
sum_lliks <- function(lliks){
  lliks[is.na(lliks)] <- -Inf
  out <- max(lliks, na.rm = TRUE) + log(sum(exp(lliks - max(lliks))))
  return(out)
}

# get weighted tip values from a stochastic map reults
get_tip_values <- function(model){
  all_tip_states <- lapply(model$simmaps, get_tip_states)
  all_joint_liks <- model$all_disc_liks + model$all_cont_liks
  weights <- exp(all_joint_liks - max(all_joint_liks))/sum(exp(all_joint_liks - max(all_joint_liks)))
  rates <- rowSums(model$solution.disc, na.rm = TRUE)
  continuous_solution <- model$solution.cont
  continuous_solution[is.na(continuous_solution)] <- 0
  parameter_table_list <- lapply(all_tip_states, function(x) index_paramers_from_tip_states(x, rates, continuous_solution))
  for(i in 1:length(parameter_table_list)){
    parameter_table_list[[i]] <- parameter_table_list[[i]] * weights[i]
  }
  tip_value_table <- Reduce("+", parameter_table_list)
  expected_values <- getExpectedValues(model, FALSE)
  expected_values <- expected_values[match(rownames(tip_value_table), expected_values$sp),]
  tip_value_table <- cbind(tip_value_table, expected_values[,c(2,3)])
  return(tip_value_table)
}

get_tip_states <- function(simmap){
  nTip <- length(simmap$tip.label)
  tip_states <- as.numeric(unlist(lapply(simmap$maps[simmap$edge[,2] <= nTip], function(x) names(x)[length(x)])))
  names(tip_states) <- simmap$tip.label[simmap$edge[,2][simmap$edge[,2] <= nTip]]
  return(tip_states)
}

index_paramers_from_tip_states <- function(tip_states, rates, continuous_solution){
  parameter_df <- data.frame(rates = rates[tip_states],
                             alpha = continuous_solution[1,tip_states],
                             sigma.sq = continuous_solution[2,tip_states],
                             theta = continuous_solution[3,tip_states],
                             row.names = names(tip_states))
  return(parameter_df)
}

# hOUwie.twostep <- function(phy, data, rate.cat, discrete_model, continuous_model, nSim=1000, root.p="yang", dual = FALSE, collapse = TRUE, root.station=FALSE, get.root.theta=FALSE, tip.fog = "none", lb_discrete_model=NULL, ub_discrete_model=NULL, lb_continuous_model=NULL, ub_continuous_model=NULL, recon=FALSE, nodes="internal", p=NULL, ip="fast", optimizer="nlopt_ln", opts=NULL, quiet=FALSE, sample_tips=FALSE, sample_nodes=TRUE, adaptive_sampling=TRUE, n_starts = 1, ncores = 1){
#   start_time <- Sys.time()
#   # if the data has negative values, shift it right - we will shift it back later
#   negative_values <- FALSE
#   if(tip.fog == "none"){
#     if(any(data[,dim(data)[2]] < 0)){
#       cat("Negative values detected... adding 50 to the trait mean for optimization purposes\n")
#       negative_values <- TRUE
#       data[,dim(data)[2]] <- data[,dim(data)[2]] + 50
#     }
#   }else{
#     if(any(data[,dim(data)[2]-1] < 0)){
#       cat("Negative values detected... adding 50 to the trait mean for optimization purposes\n")
#       negative_values <- TRUE
#       data[,dim(data)[2]-1] <- data[,dim(data)[2]-1] + 50
#     }
#   }
#   # check that tips and data match
#   # check for invariance of tip states and not that non-invariance isn't just ambiguity
#   if(!is.null(phy$node.label)){
#     if(!quiet){
#       cat("Your phylogeny had node labels, these have been removed.\n")
#     }
#     phy$node.label <- NULL
#   }
#   
#   if(ncores > n_starts){
#     cat("You have specified more cores are to be used than the number of starts. Setting ncores to be equal to the number of optimizations.\n")
#     ncores <- n_starts
#   }
#   
#   # organize the data
#   phy <- reorder.phylo(phy, "pruningwise")
#   hOUwie.dat <- organizeHOUwieDat(data, tip.fog, collapse)
#   nStates <- as.numeric(max(hOUwie.dat$data.cor[,2]))
#   nCol <- dim(data)[2] - ifelse(tip.fog == "none", 2, 3)
#   Tmax <- max(branching.times(phy))
#   tip.paths <- lapply(1:Ntip(phy), function(x) OUwie:::getPathToRoot(phy, x))
#   
#   if(class(discrete_model)[1] == "character"){
#     index.disc <- getDiscreteModel(hOUwie.dat$data.cor, discrete_model, rate.cat, dual, collapse)
#     index.disc[index.disc == 0] <- NA
#   }else{
#     index.disc <- discrete_model
#     index.disc[index.disc == 0] <- NA
#   }
#   if(class(continuous_model)[1] == "character"){
#     index.cont <- getOUParamStructure(continuous_model, "three.point", root.station, get.root.theta, nStates * rate.cat)
#   }else{
#     continuous_model[continuous_model == 0] <- NA
#     index.cont <- continuous_model
#   }
#   if(dim(index.disc)[2] > dim(index.cont)[2]){
#     stop("Not all of your discrete states have OU parameters associated with them. Please check that your discrete index matrix matches your continuous index matrix.")
#   }
#   if(dim(index.cont)[2] > dim(index.disc)[2]){
#     stop("You have specified more OU parameters than there are states in the discrete process. Please check that your discrete index matrix matches your continuous index matrix.")
#   }
#   if(class(root.p[1]) != "character"){
#     if(dim(index.disc)[2] != length(root.p)){
#       stop("You have entered a custom root prior whose length does not equal the number of states in your discrete model.")
#     }
#   }
#   
#   if(is.null(lb_continuous_model)){
#     # the lower limit of alpha is defined as a halflife of 10000% of the max tree height
#     # the lower limit of sigma is defined 10 times less than alpha
#     # the lower limit of optim is defined 10 times lower than the minimum observation
#     if(any(is.na(lb_continuous_model[1,]))){
#       lb.alpha = 1e-10
#     }else{
#       lb.alpha = 1e-10
#     }
#     lb.sigma = 1e-10
#     lb.optim = min(data[, 1+nCol+1])/10 
#     lb_continuous_model=c(lb.alpha,lb.sigma,lb.optim)
#   }
#   if(is.null(ub_continuous_model)){
#     # the upper limit of alpha is defined as a halflife of 1% of the max tree height
#     # the upper limit of sigma is defined 10 times more than alpha
#     # the upper limit of optim is defined 10 times more than the maximum observation
#     ub.alpha = log(2)/(0.01 * Tmax)
#     ub.sigma = ub.alpha
#     ub.optim = max(data[, 1+nCol+1])*10 
#     ub_continuous_model=c(ub.alpha,ub.sigma,ub.optim)
#   }
#   if(is.null(lb_discrete_model)){
#     # the minimum dwell time is defined as 100 times the max tree height
#     lb_discrete_model = 1/(Tmax*100)
#   }
#   if(is.null(ub_discrete_model)){
#     ub_discrete_model = 1/(Tmax*0.01)
#   }
#   #Ensures that weird root state probabilities that do not sum to 1 are input:
#   if(!is.null(root.p)){
#     if(!is.character(root.p)){
#       root.p <- root.p/sum(root.p)
#     }
#   }
#   time_slice <- Tmax+1
#   # the number of parameters for each process
#   n_p_trans <- max(index.disc, na.rm = TRUE)
#   n_p_alpha <- length(unique(na.omit(index.cont[1,])))
#   n_p_sigma <- length(unique(na.omit(index.cont[2,])))
#   n_p_theta <- length(unique(na.omit(index.cont[3,])))
#   n_p <- n_p_trans + n_p_alpha + n_p_sigma + n_p_theta
#   
#   # an internal data structure (internodes liks matrix) for the dev function
#   edge_liks_list <- getEdgeLiks(phy, hOUwie.dat$data.cor, nStates, rate.cat, time_slice)
#   
#   # default MLE search options
#   if(is.null(opts)){
#     if(optimizer == "nlopt_ln"){
#       opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5)
#     }
#     if(optimizer == "nlopt_gn"){
#       opts <- list("algorithm"="NLOPT_GN_DIRECT_L", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5)
#     }
#     if(optimizer == "sann"){
#       opts <- list(max.call=1000, smooth=FALSE)
#     }
#   }
#   # a global matrix to contain likelihoods so that identical parameters return identical likelihoods
#   if(is.null(opts$maxeval) | is.null(opts$max.call)){
#     max.its <- 1000
#   }else{
#     max.its <- as.numeric(opts$maxeval)
#   }
#   setDTthreads(threads=1)
#   tmp.df <- data.frame(matrix(c(0, rep(1e5, n_p)), byrow = TRUE, ncol = n_p+1, nrow = max.its))
#   global_liks_mat <- as.data.table(tmp.df)
#   
#   # p is organized into 2 groups with the first set being corHMM and the second set being OUwie
#   # organized as c(trans.rt, alpha, sigma.sq, theta)
#   # evaluate likelihood
#   if(!is.null(p)){
#     if(!quiet){
#       cat("Calculating likelihood from a set of fixed parameters.\n")
#       print(p)
#     }
#     if(max(index.cont, na.rm = TRUE) + max(index.disc, na.rm = TRUE) != length(p)){
#       message <- paste0("The number of parameters does not match the number required by the model structure. You have supplied ", length(p), ", but the model structure requires ", max(index.cont, na.rm = TRUE) + max(index.disc, na.rm = TRUE), ".")
#       stop(message, call. = FALSE)
#     }
#     out<-NULL
#     pars <- out$solution <- log(p)
#   }else{
#     out<-NULL
#     lower = log(c(rep(lb_discrete_model, n_p_trans), 
#                   rep(lb_continuous_model[1], length(unique(na.omit(index.cont[1,])))), 
#                   rep(lb_continuous_model[2], length(unique(na.omit(index.cont[2,])))), 
#                   rep(lb_continuous_model[3], length(unique(na.omit(index.cont[3,]))))))
#     upper = log(c(rep(ub_discrete_model, n_p_trans), 
#                   rep(ub_continuous_model[1], length(unique(na.omit(index.cont[1,])))), 
#                   rep(ub_continuous_model[2], length(unique(na.omit(index.cont[2,])))), 
#                   rep(ub_continuous_model[3], length(unique(na.omit(index.cont[3,]))))))
#     # cat(c("TotalLnLik", "DiscLnLik", "ContLnLik"), "\n")
#     # check for user input initial parameters 
#     if(is.character(ip)){
#       if(rate.cat > 1){
#         bin_index <- cut(hOUwie.dat$data.ou[,3], rate.cat, labels = FALSE)
#         combos <- expand.grid(1:max(hOUwie.dat$data.cor[,2]), 1:rate.cat)
#         disc_tips <- vector("numeric", length(phy$tip.label))
#         for(i in 1:dim(combos)[1]){
#           disc_tips[hOUwie.dat$data.cor[,2] == combos[i,1] & bin_index == combos[i,2]] <- i
#         }
#       }else{
#         disc_tips <- hOUwie.dat$data.cor[,2]
#       }
#       starts.alpha <- rep(log(2)/Tmax, n_p_alpha)
#       # starts.sigma <- rep(var(hOUwie.dat$data.ou[,3]), n_p_sigma)
#       starts.sigma <- rep(log(2)/Tmax, n_p_sigma)
#       start.theta <- getIP.theta(hOUwie.dat$data.ou[,3], disc_tips, index.cont[3,])
#       start.cor <- rep(10/sum(phy$edge.length), n_p_trans)
#       starts.basic = c(start.cor, starts.alpha, starts.sigma, start.theta)
#       cat("\nFitting the discrete model to discrete data...\n")
#       discrete_fit <- corHMM(phy=phy, data=hOUwie.dat$data.cor, rate.cat=rate.cat, rate.mat=index.disc, node.states="none", opts = opts)
#       cat("\nGenerating", nSim, "simmaps and optimizing the continuous model to each.")
#       simmap_list <- makeSimmap(tree=phy, data=hOUwie.dat$data.cor, model=discrete_fit$solution, rate.cat=1, nSim=nSim, nCores=1)
#       continuous_fit <- mclapply(simmap_list, function(x) nloptr(x0 = log(c(starts.alpha, starts.sigma, start.theta)), eval_f = OUwie.basic.dev, lb=lower[-seq(n_p_trans)], ub=upper[-seq(n_p_trans)], opts=opts, phy = x, data = hOUwie.dat$data.ou, tip.fog = tip.fog, index.cont = index.cont, tip.paths = tip.paths), mc.cores=ncores)
#     }
#   }
#   return(list(discrete_fit, continuous_fit))
# }

