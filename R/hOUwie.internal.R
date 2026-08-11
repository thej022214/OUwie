##### Main internal functions ##### 
withHOUwieSeed <- function(seed, code){
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if(had_seed){
    previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if(had_seed){
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    }else if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

# (alpha, sigma^2) <-> (half-life, stationary variance) is a bijection, and on the log
# scale the optimizer actually works on it is a shear: log sigma^2 = log V + log 2 +
# log alpha. The OU likelihood ridges along alpha while sigma^2/(2 alpha) is
# comparatively well determined (Cressler, Butler & King 2015, Fig. 7), so shearing
# puts the identified and the unidentified direction on separate axes rather than
# smearing both across both.
#
# Which alpha governs each sigma^2 is what decides whether the reparameterized model is
# still the model asked for. Since sigma^2_j = V_j * 2 alpha_j, a sigma^2 shared across
# regimes stays shared only if alpha is constant across those same regimes, i.e. the
# sigma^2 partition refines the alpha partition. OUM, OUMV and OUMVA all satisfy that.
# OUMA does not: alpha free per regime with one shared sigma^2 turns into the constraint
# V_j * alpha_j = constant, which no index pattern can express. Returns NULL when the
# reparameterization is unavailable, including when there is no alpha at all (BM), where
# there is no stationary variance to speak of.
alphaIndexForSigma <- function(index.cont){
  if(!length(stats::na.omit(index.cont[1,]))){
    return(NULL)
  }
  sigma_ids <- sort(unique(stats::na.omit(index.cont[2,])))
  if(!length(sigma_ids)){
    return(NULL)
  }
  governing <- vapply(sigma_ids, function(id){
    alphas <- unique(index.cont[1, which(index.cont[2,] == id)])
    if(length(alphas) != 1L || is.na(alphas)) NA_integer_ else as.integer(alphas)
  }, integer(1))
  if(anyNA(governing)){
    return(NULL)
  }
  stats::setNames(governing, sigma_ids)
}

# index.cont numbers the alpha, sigma^2 and theta parameters in one sequence, and the
# fitted vector holds them in that order after the transition rates, so a parameter's
# index in index.cont is its position in pars once shifted past the discrete block.
stationaryParPositions <- function(n_p_trans, index.cont){
  governing <- alphaIndexForSigma(index.cont)
  list(sigma = n_p_trans + as.integer(names(governing)),
       alpha_of_sigma = n_p_trans + as.integer(governing),
       alpha = n_p_trans + sort(unique(as.integer(stats::na.omit(index.cont[1,])))))
}

toStationaryPars <- function(pars, n_p_trans, index.cont){
  at <- stationaryParPositions(n_p_trans, index.cont)
  # sigma^2 is divided by the alpha it is paired with before alpha itself is replaced
  pars[at$sigma] <- pars[at$sigma]/(2 * pars[at$alpha_of_sigma])
  pars[at$alpha] <- log(2)/pars[at$alpha]
  pars
}

fromStationaryPars <- function(pars, n_p_trans, index.cont){
  at <- stationaryParPositions(n_p_trans, index.cont)
  # and alpha is recovered first here, because the stationary variances need it
  pars[at$alpha] <- log(2)/pars[at$alpha]
  pars[at$sigma] <- pars[at$sigma] * 2 * pars[at$alpha_of_sigma]
  pars
}

# Categorical draw by inverse CDF, returning 0L when no state has positive weight so
# the caller can bail rather than trap an error out of sample(). getInternodeStateSample
# inlines this rather than calling it: it draws once per internode of every edge of
# every history, where even the call overhead is worth removing.
drawFromWeights <- function(weights){
  cumulative <- cumsum(weights)
  total <- cumulative[length(cumulative)]
  if(!isTRUE(total > 0)){
    return(0L)
  }
  sum(stats::runif(1L) * total > cumulative) + 1L
}

makeCommonRandomObjective <- function(objective, seed){
  force(objective)
  force(seed)
  function(p){
    withHOUwieSeed(seed, objective(p))
  }
}

# Each start optimises against its own stream of sampled histories, so the objective
# a start finishes on is the value of that start's particular draw. Ranking starts on
# those values maximises over the random draw as well as over the parameters, and the
# winner's optimism grows with the number of starts. Re-scoring every solution under
# one seed shared by all of them separates "which start found the best parameters"
# from "which start drew the luckiest histories". Returns NULL when nothing could be
# scored, so the caller can fall back to its own ranking.
selectStartOnCommonSeed <- function(solutions, objective, seed){
  scores <- vapply(solutions, function(p){
    value <- try(withHOUwieSeed(seed, objective(p)), silent = TRUE)
    # 1e10 is the penalty hOUwie.dev returns for a numerical failure, so a solution
    # that fails under the shared seed is unrankable rather than merely bad.
    if(inherits(value, what="try-error") || !is.finite(value) || value >= 1e10){
      NA_real_
    }else{
      as.numeric(value)
    }
  }, numeric(1))
  if(all(is.na(scores))){
    return(NULL)
  }
  list(index = which.min(scores), scores = scores)
}

# A sampled history is a list of state sequences, one per edge. Every draw on a given
# tree is cut at the same internodes, so the segment lengths are constant across draws
# and carry no identifying information; one byte per state is then unambiguous without
# any separator, which is what makes the collisions this used to guard against
# impossible. States are small positive integers, so no byte is nul.
stateSampleId <- function(state.sample){
  if(length(state.sample) == 0L){
    return("")
  }
  states <- unlist(state.sample, use.names = FALSE)
  if(max(states) > 255L){
    return(paste(states, collapse = ","))
  }
  rawToChar(as.raw(states))
}

# Compile topology and discretisation metadata that is constant across every
# likelihood evaluation in one hOUwie fit. Paths are built in one cladewise
# traversal instead of searching the edge matrix once per node.
getHOUwieTreePlan <- function(phy, edge_liks_list, all.paths = NULL){
  nTip <- length(phy$tip.label)
  nEdge <- nrow(phy$edge)
  cladewise.order <- ape::reorder.phylo(phy, "cladewise", index.only = TRUE)

  if(is.null(all.paths)){
    all.paths <- vector("list", nTip + phy$Nnode)
    root <- nTip + 1L
    all.paths[[root]] <- integer(0)
    for(edge.i in cladewise.order){
      ancestor <- phy$edge[edge.i, 1]
      descendant <- phy$edge[edge.i, 2]
      all.paths[[descendant]] <- c(edge.i, all.paths[[ancestor]])
    }
  }

  edges.by.ancestor <- split(seq_len(nEdge), phy$edge[, 1])
  ancestor.match <- match(phy$edge[, 2], as.numeric(names(edges.by.ancestor)))
  child.edges <- lapply(ancestor.match, function(i){
    if(is.na(i)) integer(0) else edges.by.ancestor[[i]]
  })

  number.of.nodes.per.edge <- vapply(edge_liks_list, nrow, integer(1))
  number.of.edges.per.edge <- number.of.nodes.per.edge - 1L
  reduced.edge.length <- phy$edge.length / number.of.edges.per.edge

  list(
    all.paths = all.paths,
    tip.paths = all.paths[seq_len(nTip)],
    cladewise.order = cladewise.order,
    rev.pruning.order = rev(ape::reorder.phylo(phy, "pruningwise",
                                               index.only = TRUE)),
    root.edges = which(phy$edge[, 1] == nTip + 1L),
    child.edges = child.edges,
    edge.index = phy$edge,
    number.of.nodes.per.edge = number.of.nodes.per.edge,
    number.of.edges.per.edge = number.of.edges.per.edge,
    reduced.edge.length = reduced.edge.length,
    map.template = mapply(function(length, segments)
                            rep(length, segments),
                          length = reduced.edge.length / 2,
                          segments = number.of.edges.per.edge * 2,
                          SIMPLIFY = FALSE),
    basic.edges = cbind(seq_len(nEdge), phy$edge, 0, 0)
  )
}

# hOUwie.dev draws its own maps, so a draw that leaves no usable map is recoverable: the
# optimizer only needs a finite penalty to be able to move away from the point. split.liks
# is the reporting call, and its result is handed to getHouwieObj, which indexes it as a
# list, so there is no scalar it can be given and no likelihood to report.
houwieDevFailure <- function(split.liks, reason){
  if(split.liks){
    stop(paste0("The hOUwie likelihood cannot be reported at these parameters because ", reason, ". Check root.p and the discrete model against the states present in your data."), call. = FALSE)
  }
  return(1e10)
}

hOUwie.dev <- function(p, phy, data, rate.cat, tip.fog,
                       index.disc, index.cont, root.p,
                       edge_liks_list, nSim, all.paths=NULL, 
                       sample_nodes=TRUE, split.liks=FALSE, 
                       global_liks_mat=NULL, diagn_msg=FALSE,
                       tree.plan=NULL, algorithm="sampling",
                       resolution=1L, max_components=Inf, tolerance=0,
                       history="midpoint", proposal_temper=1,
                       proposal_defensive=0.5,
                       parameterization="alpha-sigma"){
  if(is.null(tree.plan)){
    tree.plan <- getHOUwieTreePlan(phy, edge_liks_list, all.paths)
  }
  all.paths <- tree.plan$all.paths
  tip.paths <- tree.plan$tip.paths
  p <- exp(p)
  # check if these parameters exist in the global matrix
  # set(global_liks_mat, i = as.integer(1),  j = 1:4, value=as.list(c(0, p)))
  if(!is.null(global_liks_mat)){
    # a row matches only when every parameter matches. summing the signed differences
    # instead would call two different parameter vectors identical whenever their
    # deviations happen to cancel, and hand back the wrong cached likelihood.
    liks_match_vector <- colSums(t(global_liks_mat[,-1]) != p) == 0
    match_row <- which(liks_match_vector)[1]
    cached_llik_houwie <- if(is.na(match_row)) NA_real_ else as.numeric(global_liks_mat[match_row, 1])
    if(!split.liks){
      if(!is.na(match_row)){
        # print(cached_llik_houwie)
        # print(p)
        return(-cached_llik_houwie)
      }
    }
  }

  k <- max(index.disc, na.rm = TRUE)
  p.mk <- p[1:k]
  p.ou <- p[(k+1):length(p)]
  Rate.mat <- matrix(1, 3, dim(index.disc)[2])
  index.cont[is.na(index.cont)] <- max(index.cont, na.rm = TRUE) + 1
  Rate.mat[] <- c(p.ou, 1e-10)[index.cont]
  alpha = Rate.mat[1,]
  sigma.sq = Rate.mat[2,]
  theta = Rate.mat[3,]
  if(parameterization == "halflife-variance"){
    # the alpha and sigma rows carry the half-life and the stationary variance on this
    # path. alpha is shared across regimes, so the per-regime stationary variances map
    # back elementwise. See toStationaryPars().
    alpha <- log(2)/alpha
    sigma.sq <- sigma.sq * 2 * alpha
  }
  rate <- index.disc
  rate[is.na(rate)] <- k + 1
  Q <- matrix(0, dim(rate)[1], dim(rate)[2])
  Q[] <- c(p.mk, 0)[rate]
  diag(Q) <- -rowSums(Q)
  edge_liks_list_init <- edge_liks_list
  if(algorithm == "pruning"){
    # The pruning path sums over regime histories exactly rather than sampling
    # them, so none of the proposal machinery below applies. Tip constraints are
    # read off edge_liks_list so ambiguous coding and hidden states behave the
    # same way they do on the sampling path.
    tip_values <- data[match(phy$tip.label, data[,1]), 3]
    llik_pruned <- try(houwiePruningLik(
      phy = phy,
      tip_state_sets = getTipStateSets(phy, edge_liks_list_init,
                                       dim(index.disc)[2]),
      tip_values = tip_values, Q = Q, alpha = alpha, sigma.sq = sigma.sq,
      theta = theta, root.p = root.p, resolution = resolution,
      max_components = max_components, tolerance = tolerance), silent = TRUE)
    if(inherits(llik_pruned, what="try-error") || !is.finite(llik_pruned)){
      return(houwieDevFailure(split.liks, "the pruned likelihood could not be computed"))
    }
    llik_pruned <- as.numeric(llik_pruned)
    if(split.liks){
      # the discrete marginal is available exactly from the same Q, so it is
      # reported rather than approximated; the continuous part is not separable
      # once the history has been integrated out, so it is not invented
      return(list(TotalLik = llik_pruned, DiscLik = NA_real_,
                  ContLik = NA_real_, llik_discrete = numeric(0),
                  llik_continuous = numeric(0), simmaps = NULL,
                  unsorted_lliks_df = NULL))
    }
    if(!is.null(global_liks_mat)){
      new_row <- which(global_liks_mat$X1 == 0)[1]
      if(!is.na(new_row)){
        set(global_liks_mat, as.integer(new_row), names(global_liks_mat),
            as.list(c(llik_pruned, p)))
      }
    }
    if(diagn_msg){
      print(c(round(llik_pruned, 2), round(p, 2)))
    }
    return(-llik_pruned)
  }
  # The continuous factors are built from the messages below each node and are
  # applied after the discrete conditionals have been propagated, so they are
  # computed here and used further down rather than folded into edge_liks_list.
  continuous_node_factors <- NULL
  if(sample_nodes){
    tip_values <- data[match(phy$tip.label, data[,1]), 3]
    continuous_node_factors <- try(getContinuousNodeFactors(
      phy, tip_values, Q, alpha, sigma.sq, theta, temper = proposal_temper), silent = TRUE)
    if(inherits(continuous_node_factors, what="try-error")){
      return(houwieDevFailure(split.liks, "the continuous messages used to shape the proposal could not be computed"))
    }
  }
  # get the condtional probabilities based on the discrete values
  for(recon_index in 1:length(edge_liks_list)){
    edge_liks_list[[recon_index]] <- edge_liks_list[[recon_index]] * edge_liks_list_init[[recon_index]]
  }
  conditional_probs <- getConditionalInternodeLik(phy, Q, edge_liks_list)
  defensive_probs <- NULL
  if(!is.null(continuous_node_factors)){
    # the untilted conditionals are kept and become the defensive component of the
    # proposal mixture, which is what stops a saturated tilt from collapsing the draw
    defensive_probs <- conditional_probs
    conditional_probs <- applyContinuousNodeFactors(conditional_probs, phy,
                                                    continuous_node_factors)
  }
  root_liks <- getRootLiks(conditional_probs, Q, root.p)
  if(is.null(root_liks)){
    return(houwieDevFailure(split.liks, "no state at the root has a defined conditional probability"))
  }
  # initial sample
  # sample mappings based on the conditional probabilites (also calculating some time saving probabilities from transitions to and from particular states)
  # Bridge histories are drawn per parameter evaluation, so the uniformization
  # tables are built here and shared by every branch of every draw.
  uni <- NULL
  if(history == "bridge"){
    uni <- makeUniformizationTable(Q, max(tree.plan$reduced.edge.length))
    if(is.null(uni)){
      return(houwieDevFailure(split.liks, "the transition rates are too large for an exact regime path to be drawn"))
    }
  }
  internode_maps_and_discrete_probs <- getInternodeMap(phy, Q, conditional_probs$edge_liks_list, conditional_probs$root_state, root_liks, nSim, check_vector = NA, max.attempts=nSim*2, tree.plan=tree.plan, unique_only=FALSE, uni=uni,
    defensive_edge_liks_list = if(is.null(defensive_probs)) NULL else defensive_probs$edge_liks_list,
    defensive_root_state = if(is.null(defensive_probs)) NULL else defensive_probs$root_state,
    defensive_lambda = proposal_defensive)
  internode_maps <- internode_maps_and_discrete_probs$maps
  internode_samples <- internode_maps_and_discrete_probs$state_samples
  log_proposal <- internode_maps_and_discrete_probs$log_proposal
  mapping_ids <- internode_maps_and_discrete_probs$mapping_ids
  check_vector <- mapping_ids
  # The draw is now taken with replacement, so it always reaches nSim unless a
  # history failed outright. The Q * 100 top-up that used to run here proposed
  # from a different distribution whose density was never recorded, which the
  # importance weights cannot absorb.
  # calculte the discrete probabilities based on the given Q matrix (Pij already calculated)
  discrete_probs <- lapply(internode_samples, function(x) getStateSampleProb(state_sample = x, Pij = internode_maps_and_discrete_probs$Pij, root_liks = root_liks, root_edges = internode_maps_and_discrete_probs$root_edges))
  llik_discrete <- unlist(discrete_probs)
  failed_maps <- discrete_probs == -Inf
  llik_discrete <- llik_discrete[!failed_maps]
  log_proposal <- log_proposal[!failed_maps]
  mapping_ids <- mapping_ids[!failed_maps]
  if(length(llik_discrete) == 0){
    return(houwieDevFailure(split.liks, "every sampled map has probability zero under the discrete model"))
  }
  # Keep histories as compact edge-segment lists while evaluating the objective.
  # Full simmap trees and mapped.edge matrices are only needed by downstream output.
  compact.maps <- internode_maps[!failed_maps]
  # Left to itself OUwie.basic recovers the regime set with unique(unlist(lapply(map,
  # names))), which walks every segment of every edge in character space once per
  # draw. The regimes are known here and are the same for every history, so they are
  # handed over instead. Passing all of them rather than only the ones a particular
  # history visits also skips the subsetting of Rate.mat and pars, and leaves regime
  # numbers indexing those directly instead of by position in the visited subset.
  map.states <- as.character(seq_along(theta))
  continuous.lik <- function(current.map){
    OUwie.basic(phy, data, simmap.tree=TRUE, scaleHeight=FALSE,
                alpha=alpha, sigma.sq=sigma.sq, theta=theta,
                algorithm="three.point", tip.paths=tip.paths, tip.fog=tip.fog,
                map=current.map, map.states=map.states, tree.plan=tree.plan)
  }
  simmaps <- NULL
  # if there is no character dependence the map has no influence on continuous likleihood
  character_dependence_check <- all(apply(index.cont, 1, function(x) length(unique(x)) == 1))
  if(character_dependence_check){
    llik_continuous <- continuous.lik(compact.maps[[1]])
    llik_continuous <- rep(llik_continuous, length(compact.maps))
  }else{
    # repeated histories are evaluated once and reused, which is the whole of the
    # saving that dropping them used to buy, without changing the draw. bridge
    # histories share only their node states, so two draws with the same id have
    # different paths and different continuous likelihoods, and reusing one for
    # the other would be wrong
    if(is.null(uni)){
      unique_index <- match(mapping_ids, unique(mapping_ids))
      unique_liks <- unlist(lapply(compact.maps[!duplicated(unique_index)],
                                   continuous.lik))
      llik_continuous <- unique_liks[unique_index]
    }else{
      llik_continuous <- unlist(lapply(compact.maps, continuous.lik))
    }
  }
  # combine probabilities being careful to avoid underflow
  llik_houwies <- llik_discrete + llik_continuous
  # Importance sampling estimate of L = sum_h P(h | Q) f(x | h). The histories
  # were drawn from q, not from P(. | Q), so each term carries the weight
  # P(h | Q) / q(h); dividing by the number of draws is what makes this an
  # average rather than the subset sum it used to be. Summing the terms without
  # the weight and without the 1/n effectively weights each history by its own
  # proposal probability, which rewards parameter values that concentrate the
  # sampling distribution on few histories.
  log_terms <- llik_discrete - log_proposal + llik_continuous
  n_draws <- length(log_terms)
  finite_terms <- log_terms[is.finite(log_terms)]
  if(!length(finite_terms)){
    return(houwieDevFailure(split.liks, "no sampled map has a finite weight"))
  }
  max_term <- max(finite_terms)
  llik_houwie <- max_term + log(sum(exp(log_terms - max_term))) - log(n_draws)
  # Effective sample size of the normalized weights. This is the number that says
  # whether the estimate can be believed: at n_eff of 1 a single history carries
  # the whole likelihood and the proposal is doing nothing.
  normalized <- exp(log_terms - max_term)
  normalized[!is.finite(normalized)] <- 0
  ess <- sum(normalized)^2 / sum(normalized^2)
  # the weight formula returns exactly n when every draw is identical, so a collapsed
  # proposal reports perfect efficiency. A draw of k distinct histories carries at most
  # k effective samples whatever the weights say. Bridge histories sharing node states
  # are still different paths, so the cap only applies where an id really is the history.
  if(is.null(uni)){
    ess <- min(ess, length(unique(mapping_ids)))
  }

  unsorted_lliks_df <- data.frame(llik_discrete=llik_discrete, llik_continuous=llik_continuous)
  # the maps themselves are still reported best-first, which only affects output
  sorted_likelihoods <- sort(llik_houwies, decreasing = TRUE, index.return = TRUE)
  selected.map.index <- c(na.omit(sorted_likelihoods$ix[1:nSim]))
  compact.maps <- compact.maps[selected.map.index]
  llik_discrete_summed <- max(llik_discrete) + log(sum(exp(llik_discrete - max(llik_discrete))))
  llik_continuous_summed <- max(llik_continuous) + log(sum(exp(llik_continuous - max(llik_continuous))))
  if(split.liks){
    if(is.null(simmaps)){
      simmaps <- getMapFromSubstHistory(compact.maps, phy)
    }
    # expected_vals <- lapply(simmaps, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, tip.fog=tip.fog,return.expected.vals=TRUE))
    # expected_vals <- colSums(do.call(rbind, expected_vals) * exp(llik_houwies - max(llik_houwies))/sum(exp(llik_houwies - max(llik_houwies))))
    # report the value the optimizer actually saw at this point rather than a fresh
    # draw over a new set of maps, which would not generally reproduce it
    if(!is.null(global_liks_mat) && !is.na(cached_llik_houwie)){
      llik_houwie <- cached_llik_houwie
    }
    return(list(TotalLik = llik_houwie, DiscLik = llik_discrete_summed, ContLik = llik_continuous_summed, llik_discrete=llik_discrete, llik_continuous=llik_continuous, simmaps=simmaps, unsorted_lliks_df=unsorted_lliks_df, ess=ess))
  }
  # every map underflowing leaves max() + log(sum(exp())) as NaN, which the optimizer
  # cannot order against anything. a failed evaluation has to leave here as the same
  # finite penalty every other failure returns. the sentinel is a scalar, so it can only
  # be returned on the optimizer's path - split.liks callers hand the result to
  # getHouwieObj, which indexes it as a list
  if(!is.finite(llik_houwie)){
    return(1e10)
  }
  if(!is.null(global_liks_mat)){
    # no free row means the optimizer ran past the cap the table was sized for; the
    # cache just stops growing rather than silently discarding the write
    new_row <- which(global_liks_mat$X1 == 0)[1]
    if(!is.na(new_row)){
      set(global_liks_mat, as.integer(new_row), names(global_liks_mat), as.list(c(llik_houwie, p)))
    }
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
                             sample_nodes=TRUE, split.liks=FALSE,
                             global_liks_mat=NULL, diagn_msg=FALSE){
  tip.paths <- all.paths[1:length(simmaps[[1]]$tip.label)]
  p <- exp(p)
  # check if these parameters exist in the global matrix
  # set(global_liks_mat, i = as.integer(1),  j = 1:4, value=as.list(c(0, p)))
  if(!is.null(global_liks_mat)){
    # a row matches only when every parameter matches. summing the signed differences
    # instead would call two different parameter vectors identical whenever their
    # deviations happen to cancel, and hand back the wrong cached likelihood.
    liks_match_vector <- colSums(t(global_liks_mat[,-1]) != p) == 0
    match_row <- which(liks_match_vector)[1]
    cached_llik_houwie <- if(is.na(match_row)) NA_real_ else as.numeric(global_liks_mat[match_row, 1])
    if(!split.liks){
      if(!is.na(match_row)){
        # print(cached_llik_houwie)
        # print(p)
        return(-cached_llik_houwie)
      }
    }
  }

  k <- max(index.disc, na.rm = TRUE)
  p.mk <- p[1:k]
  p.ou <- p[(k+1):length(p)] 
  Rate.mat <- matrix(1, 3, dim(index.disc)[2])
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
  # the maps are supplied by the caller and are never resampled, so every one of them
  # having probability zero is a property of the model rather than of the point being
  # evaluated. no parameter value recovers it, and neither caller can be given a number:
  # split.liks hands the result to getHouwieObj, which indexes it as a list, and the
  # optimizer would search a surface that is flat everywhere.
  if(length(llik_discrete) == 0){
    stop("Every supplied stochastic map has probability zero under this model, so the hOUwie likelihood is undefined. Either root.p gives zero probability to the state the maps have at the root, or the discrete model disallows a transition that the maps contain.", call. = FALSE)
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
    # report the value the optimizer actually saw at this point rather than a fresh
    # draw over a new set of maps, which would not generally reproduce it
    if(!is.null(global_liks_mat) && !is.na(cached_llik_houwie)){
      llik_houwie <- cached_llik_houwie
    }
    return(list(TotalLik = llik_houwie, DiscLik = llik_discrete_summed, ContLik = llik_continuous_summed, llik_discrete=llik_discrete, llik_continuous=llik_continuous, simmaps=simmaps))
  }
  # every map underflowing leaves max() + log(sum(exp())) as NaN, which the optimizer
  # cannot order against anything. a failed evaluation has to leave here as the same
  # finite penalty every other failure returns. the sentinel is a scalar, so it can only
  # be returned on the optimizer's path - split.liks callers hand the result to
  # getHouwieObj, which indexes it as a list
  if(!is.finite(llik_houwie)){
    return(1e10)
  }
  if(!is.null(global_liks_mat)){
    # no free row means the optimizer ran past the cap the table was sized for; the
    # cache just stops growing rather than silently discarding the write
    new_row <- which(global_liks_mat$X1 == 0)[1]
    if(!is.na(new_row)){
      set(global_liks_mat, as.integer(new_row), names(global_liks_mat), as.list(c(llik_houwie, p)))
    }
  }
  if(diagn_msg){
    print(c(round(llik_houwie, 2), round(llik_discrete_summed, 2), round(llik_continuous_summed, 2), round(p, 2)))
  }
  # print(p)
  return(-llik_houwie)
}

# internal for houwie.thorough
# model average the fitted parameter matrices of a set of models so that they can be
# used as shared starting values. getModelAvgParams is reserved for tip averaging.
getModelAvgStartingPars <- function(model.list, type="BIC", BM_alpha_treatment="zero"){
  if(length(model.list) == 1){
    AICwts <- 1
  }else{
    mods_table <- getModelTable(model.list, type=type)
    AICwts <- mods_table[,grep("wt$", colnames(mods_table))]
  }
  solution_disc <- lapply(model.list, "[[", "solution.disc")
  solution_cont <- lapply(model.list, "[[", "solution.cont")
  # a transition disallowed by one model contributes nothing to the average, as does
  # alpha for a BM model when it is treated as zero
  disc_table <- do.call(rbind, lapply(solution_disc, c))
  disc_table[is.na(disc_table)] <- 0
  cont_table <- do.call(rbind, lapply(solution_cont, c))
  if(BM_alpha_treatment == "zero"){
    cont_table[is.na(cont_table)] <- 0
  }
  mod_avg_disc <- matrix(colSums(disc_table * AICwts), dim(solution_disc[[1]])[1], dim(solution_disc[[1]])[2],
                         dimnames = dimnames(solution_disc[[1]]))
  mod_avg_cont <- matrix(colSums(cont_table * AICwts), dim(solution_cont[[1]])[1], dim(solution_cont[[1]])[2],
                         dimnames = dimnames(solution_cont[[1]]))
  return(list(mod_avg_disc = mod_avg_disc, mod_avg_cont = mod_avg_cont))
}

runSingleThorough <- function(houwie_obj, new_maps, init_pars){
  hOUwie.dat <- houwie_obj$hOUwie.dat
  root.p <- houwie_obj$root.p
  tip.fog <- houwie_obj$tip.fog
  rate.cat <- houwie_obj$rate.cat
  index.disc <- houwie_obj$index.disc
  index.cont <- houwie_obj$index.cont
  # the free parameters of each model are numbered 1:k within their index matrix, so an
  # averaged parameter is the mean of the cells a particular model ties together
  n_p_trans <- max(index.disc, na.rm = TRUE)
  p_disc <- unlist(lapply(1:n_p_trans, function(x) mean(init_pars$mod_avg_disc[which(index.disc == x)])))
  n_p_cont <- max(index.cont, na.rm = TRUE)
  p_cont <- unlist(lapply(1:n_p_cont, function(x) mean(init_pars$mod_avg_cont[which(index.cont == x)])))
  ip <- c(p_disc, p_cont)
  # every start is searched on the log scale, so each has to be strictly positive on the
  # scale it is finally logged on. thetas are averaged on the original trait scale, and
  # hOUwie.fixed shifts a negative trait right and carries the theta block of ip along with
  # it, so a theta is judged after that same shift while the rates, alphas and sigmas - which
  # are never shifted - are judged as they stand. testing thetas before the shift would
  # discard the averaging for any negative trait; exempting them would let a theta that stays
  # non-positive through to log(), which only warns.
  n_p_theta <- length(unique(na.omit(index.cont[3,])))
  trait_column <- ifelse(tip.fog == "none", dim(hOUwie.dat$data.ou)[2], dim(hOUwie.dat$data.ou)[2] - 1)
  shifted_ip <- ip
  if(n_p_theta > 0){
    theta_index <- (length(ip) - n_p_theta + 1):length(ip)
    shifted_ip[theta_index] <- shifted_ip[theta_index] + getTraitShift(hOUwie.dat$data.ou[,trait_column])
  }
  # fall back on the default starting values if averaging produced something unusable
  if(any(!is.finite(ip)) | any(shifted_ip <= 0)){
    warning("Model averaged starting values were not usable for one of the models, default starting values were used instead.")
    ip <- NULL
  }
  res <- hOUwie.fixed(simmaps = new_maps, data = hOUwie.dat$data.ou, rate.cat = rate.cat, discrete_model = index.disc, continuous_model = index.cont, root.p = root.p, tip.fog = tip.fog, make_numeric = FALSE, ip = ip)
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


# `resolution` counts regime points inserted per edge and matches the argument of
# the same name on the pruning path: 1 is hOUwie's own history model, in which an
# edge whose endpoints differ switches at its midpoint. It is expressed per edge
# rather than per unit time because what the painting has to resolve is where on
# a branch a shift sits, and time_slice only ever refines branches longer than
# itself - on a tree whose branches are mostly shorter it changes nothing.
getEdgeLiks <- function(phy, data, n.traits, rate.cat, time_slice, resolution = 1L){
  edge_liks_list <- vector("list", dim(phy$edge)[1])
  nTip <- length(phy$tip.label)
  for(edge_i in 1:dim(phy$edge)[1]){
    # +2 because we slice the middle of the branch and need 2 terminal nodes (ancestor and descendent)
    n_slice <- if(resolution > 1L){
      as.integer(resolution) + 1L
    }else{
      (phy$edge.length[edge_i] %/% time_slice) + 2
    }
    edge_liks_list[[edge_i]] <- matrix(1, n_slice, n.traits * rate.cat)
    if(phy$edge[edge_i,2] <= nTip){
      tmp <- numeric(n.traits)
      species_i <- phy$tip.label[phy$edge[edge_i,2]]
      if(!species_i %in% data[,1]){
        next
      }
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


# defensive_edge_liks_list / defensive_root_state give a second proposal to mix with.
# Both endpoints of the tilt are bad on their own -- untilted is far too diffuse (ESS
# near 1 where the trait is informative), fully tilted saturates onto the MAP history --
# so histories are drawn from q = lambda*q_tilted + (1-lambda)*q_untilted and every
# history is scored under both components. The mixture density is exact, so the
# estimator stays unbiased, and the weights are bounded by 1/lambda times the better
# component's: whichever proposal would have been right, this cannot do much worse.
getInternodeMap <- function(phy, Q, edge_liks_list, root_state, root_liks, nSim,
                           check_vector=NULL, max.attempts, tree.plan=NULL,
                           unique_only=TRUE, uni=NULL,
                           defensive_edge_liks_list=NULL, defensive_root_state=NULL,
                           defensive_lambda=0.5){
  # set-up
  current.attempts <- 0
  nStates <- dim(Q)[1]
  nTip <- length(phy$tip.label)
  # a potential speedup is to calculate all Pij (bollback eq.3) for all branches first
  Pij <- array(0, c(dim(Q)[1], dim(Q)[2], length(phy$edge.length)))
  # reduced edge.lengths since we are including internodes
  if(is.null(tree.plan)){
    tree.plan <- getHOUwieTreePlan(phy, edge_liks_list)
  }
  number_of_nodes_per_edge <- tree.plan$number.of.nodes.per.edge
  number_of_edges_per_edge <- tree.plan$number.of.edges.per.edge
  reduced_edge_length <- tree.plan$reduced.edge.length
  for(i in 1:length(phy$edge.length)){
    Pij[,,i] <- expm(Q * reduced_edge_length[i])
  }
  # the probability of a descendent being in state j given starting in the row of the Pj matrix
  buildPj <- function(liks){
    out <- vector("list", length(phy$edge.length))
    for(i in seq_along(out)){
      out[[i]] <- array(0, c(dim(Q)[1], dim(Q)[2], number_of_nodes_per_edge[i]))
      for(j in 1:number_of_nodes_per_edge[i]){
        out[[i]][,,j] <- sweep(Pij[,,i], MARGIN = 2, liks[[i]][j,], '*')
      }
    }
    out
  }
  Pj <- buildPj(edge_liks_list)
  # the transitions are shared, so only the sweep against the conditionals differs
  defensive <- !is.null(defensive_edge_liks_list) && defensive_lambda < 1
  Pj_defensive <- if(defensive) buildPj(defensive_edge_liks_list) else NULL
  # simulate nSim substitution histories
  rev.pruning.order <- tree.plan$rev.pruning.order
  sub_histories <- vector("list", nSim)
  root_edges <- tree.plan$root.edges
  edge_index <- tree.plan$edge.index
  # the edges descending from each edge never change, so they are listed once here
  # rather than being searched for again on every edge of every simulated history
  child_edges <- tree.plan$child.edges
  Map_i <- tree.plan$map.template
  if(unique_only && !is.null(check_vector)){
    state_samples <- vector("list", nSim)
    sim_counter <- 0
    while(!(sim_counter >= nSim | current.attempts >= max.attempts)){
      state_sample <- try(getInternodeStateSample(Pj, root_state, root_edges, rev.pruning.order, edge_index, nStates, number_of_nodes_per_edge, child_edges), silent = TRUE)
      if(inherits(state_sample, "try-error")){
        current.attempts <- current.attempts + 1
      }else{
        current_mapping_id <- stateSampleId(state_sample)
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
    # which component each history comes from is drawn up front, so the count of draws
    # per component is itself random in the right way for the mixture density below
    from_tilted <- if(defensive) stats::runif(nSim) < defensive_lambda else rep(TRUE, nSim)
    state_samples <- lapply(seq_len(nSim), function(x){
      Pj_x <- if(from_tilted[x]) Pj else Pj_defensive
      root_x <- if(from_tilted[x]) root_state else defensive_root_state
      try(getInternodeStateSample(Pj_x, root_x, root_edges, rev.pruning.order, edge_index, nStates, number_of_nodes_per_edge, child_edges), silent = TRUE)
    })
    state_samples <- state_samples[unlist(lapply(state_samples, class)) != "try-error"]
  }
  mapping_ids <- vapply(state_samples, stateSampleId, character(1))
  # Dropping repeats turns the draw into a sample without replacement, which no
  # longer has the proposal density the importance weights assume. Repeats are
  # kept when unique_only is FALSE and deduplicated only for the evaluation of
  # the continuous likelihood, which is where the saving actually was.
  if(unique_only){
    state_samples <- state_samples[!duplicated(mapping_ids)]
    mapping_ids <- vapply(state_samples, stateSampleId, character(1))
  }
  state_samples <- state_samples[!mapping_ids == ""]
  if(is.null(uni)){
    maps <- lapply(state_samples, function(x) getMapFromStateSample(Map_i, x))
  }else{
    # a bridge that cannot be represented takes its history out of the draw
    # rather than being replaced by an approximate path, which would carry a
    # density the importance weights do not know about
    maps <- lapply(state_samples, function(x)
      getMapFromStateSampleBridge(x, phy$edge.length, uni))
    drawn <- !vapply(maps, is.null, logical(1))
    maps <- maps[drawn]
    state_samples <- state_samples[drawn]
  }
  # a history's density is its density under the mixture, not under the component it
  # happened to be drawn from, so each one is scored under both
  log_proposal <- if(defensive){
    vapply(state_samples, function(x){
      lp_tilted <- scoreInternodeStateSample(x, Pj, root_state, root_edges, rev.pruning.order, nStates)
      lp_plain <- scoreInternodeStateSample(x, Pj_defensive, defensive_root_state, root_edges, rev.pruning.order, nStates)
      terms <- c(log(defensive_lambda) + lp_tilted, log1p(-defensive_lambda) + lp_plain)
      terms <- terms[is.finite(terms)]
      if(!length(terms)) return(-Inf)
      max(terms) + log(sum(exp(terms - max(terms))))
    }, numeric(1))
  }else{
    vapply(state_samples, function(x) attr(x, "log_proposal"), numeric(1))
  }
  return(list(state_samples=state_samples, maps = maps, root_edges=root_edges,
              Pij = Pij, Pj = Pj,
              log_proposal = log_proposal,
              mapping_ids = vapply(state_samples, stateSampleId,
                                   character(1))))
}

getMapFromStateSample <- function(map, state_sample){
  for(edge_i in seq_along(map)){
    states_i <- state_sample[[edge_i]]
    n_i <- length(states_i)
    names(map[[edge_i]]) <- if(n_i > 2L){
      c(states_i[1L], rep(states_i[-c(1L, n_i)], each = 2L), states_i[n_i])
    }else{
      states_i
    }
  }
  return(map)
}

# A history is drawn root-tipward from the conditional probabilities in Pj, so
# the proposal is not the model. The log density of the draw is accumulated here
# and returned with it, because the importance weight P(h | Q) / q(h) needs it
# and it cannot be recovered afterwards: Pj carries the partial likelihoods that
# tilt the draw, and only the normalizer of each categorical step records how
# much it tilted. The draw is taken by inverse CDF rather than with sample(), whose
# argument validation costs more than the draw itself when it runs once per internode
# of every edge of every history; the cumulative sum it walks also hands back the
# normalizer the density needs, so nothing is summed twice.
getInternodeStateSample <- function(Pj, root_state, root_edge, rev.pruning.order, edge_index, nStates, number_of_nodes_per_edge, child_edges=NULL){
  # each map will have edges split into equal time portions. States are held as
  # integers so that indexing with them, and the as.character() that turns them into
  # map names downstream, both skip a double-to-integer round trip on every draw.
  state_samples <- lapply(number_of_nodes_per_edge, function(x) integer(x))
  # one block of uniforms for the whole history: the root draw plus one per internode
  # step, so the generator is entered once instead of once per state.
  uniforms <- runif(1L + sum(number_of_nodes_per_edge - 1L))
  draw_i <- 1L

  cumulative <- cumsum(root_state)
  total <- cumulative[length(cumulative)]
  if(!isTRUE(total > 0)){
    stop("no state at the root has positive conditional probability")
  }
  # a state whose weight is zero spans an empty interval of the cumulative sum, so it
  # can never be drawn and log(weight) below is never -Inf
  root_sample <- sum(uniforms[draw_i] * total > cumulative) + 1L
  draw_i <- draw_i + 1L
  log_proposal <- log(root_state[root_sample]) - log(total)
  for(i in root_edge){
    state_samples[[i]][1] <- root_sample
  }
  # sample the nodes along a branch the last dec node goes into the next map
  for(edge_i in rev.pruning.order){
    from <- state_samples[[edge_i]][1]
    count <- 2
    n_inter_nodes <- length(state_samples[[edge_i]])
    for(inter_edge_i in (n_inter_nodes-1):1){
      step_weights <- Pj[[edge_i]][from,,inter_edge_i]
      cumulative <- cumsum(step_weights)
      total <- cumulative[nStates]
      if(!isTRUE(total > 0)){
        stop("no state has positive conditional probability at an internode")
      }
      if(draw_i > length(uniforms)){
        uniforms <- c(uniforms, runif(length(uniforms)))
      }
      to <- sum(uniforms[draw_i] * total > cumulative) + 1L
      draw_i <- draw_i + 1L
      log_proposal <- log_proposal + log(step_weights[to]) - log(total)
      from <- state_samples[[edge_i]][count] <- to
      count <- count + 1
    }
    if(is.null(child_edges)){
      anc_edge <- which(edge_index[edge_i,2] == edge_index[,1])
    }else{
      anc_edge <- child_edges[[edge_i]]
    }
    for(i in anc_edge){
      state_samples[[i]][1] <- to
    }
  }
  attr(state_samples, "log_proposal") <- log_proposal
  return(state_samples)
}


# Density of an already drawn history under a given set of internode conditionals.
# getInternodeStateSample accumulates this while it draws, but a defensive mixture has
# to score each history under the component it did *not* come from as well, so the same
# recursion is available here with the states held fixed. The traversal order is
# irrelevant when nothing is being sampled; it is kept the same only so the two agree
# term for term.
scoreInternodeStateSample <- function(state_samples, Pj, root_state, root_edge,
                                      rev.pruning.order, nStates){
  root_sample <- state_samples[[root_edge[1]]][1]
  total <- sum(root_state)
  if(!isTRUE(total > 0) || !isTRUE(root_state[root_sample] > 0)){
    return(-Inf)
  }
  log_density <- log(root_state[root_sample]) - log(total)
  for(edge_i in rev.pruning.order){
    from <- state_samples[[edge_i]][1]
    count <- 2L
    n_inter_nodes <- length(state_samples[[edge_i]])
    for(inter_edge_i in (n_inter_nodes-1):1){
      step_weights <- Pj[[edge_i]][from,,inter_edge_i]
      to <- state_samples[[edge_i]][count]
      step_total <- sum(step_weights)
      if(!isTRUE(step_total > 0) || !isTRUE(step_weights[to] > 0)){
        return(-Inf)
      }
      log_density <- log_density + log(step_weights[to]) - log(step_total)
      from <- to
      count <- count + 1L
    }
  }
  log_density
}

# get path probability internal
# The step probabilities are read in one matrix-indexing pass. Walking the path with
# path_states <- path_states[-1] instead reallocated the path on every step, which made
# a linear traversal quadratic in the number of internodes on the edge.
getPathStateProb <- function(path_states, p_mat){
  n <- length(path_states)
  if(n < 2L){
    return(0)
  }
  sum(log(p_mat[path_states[-n] + (path_states[-1L] - 1L) * dim(p_mat)[1L]]))
}


# Every from/to step of every edge is turned into a linear position in Pij and read in
# one pass. The per-edge loop cost a matrix slice of Pij per edge per draw, which is
# where most of the discrete probability time went.
getStateSampleProb <- function(state_sample, Pij, root_liks, root_edges){
  root_sample <- state_sample[[root_edges[1]]][1] # the root sample
  n_steps <- lengths(state_sample) - 1L
  if(all(n_steps <= 0L)){
    return(log(root_liks[root_sample]))
  }
  nStates <- dim(Pij)[1L]
  from <- unlist(lapply(state_sample, function(x) x[-length(x)]), use.names = FALSE)
  to <- unlist(lapply(state_sample, function(x) x[-1L]), use.names = FALSE)
  edge <- rep.int(seq_along(state_sample), pmax(n_steps, 0L))
  positions <- from + (to - 1L) * nStates + (edge - 1L) * nStates * nStates
  llik <- sum(log(Pij[positions])) + log(root_liks[root_sample])
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
  # the row names are a property of the tree, not of the map, so they are built once
  # here rather than once per simmap
  RowNames <- paste(phy$edge[,1], phy$edge[,2], sep = ",")
  mapped.edge <- lapply(maps, function(x) convertSubHistoryToEdge(phy, x, RowNames))
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
OUwie.basic <- function(phy, data, simmap.tree=TRUE, root.age=NULL,
                        scaleHeight=FALSE, root.station=FALSE,
                        get.root.theta=FALSE, shift.point=0.5, alpha,
                        sigma.sq, theta, tip.fog="none",
                        algorithm="three.point", tip.paths=NULL,
                        return.expected.vals=FALSE, map=NULL,
                        map.states=NULL, tree.plan=NULL){
  compact.map <- !is.null(map)
  if(!compact.map){
    map <- phy$maps
  }
  # organize tip states based on what the simmap suggests
  nTip <- length(phy$tip.label)
  if(!compact.map){
    mapping <- unlist(lapply(map, function(x) names(x[length(x)])))
    TipStates <- mapping[match(match(data[,1], phy$tip.label), phy$edge[,2])]
    data[,2] <- TipStates
  }
  
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
  if(is.null(map.states)){
    if(compact.map){
      map.states <- unique(unlist(lapply(map, names), use.names = FALSE))
    }else{
      map.states <- colnames(phy$mapped.edge)
    }
    # both sources list the regimes in the order the map happens to mention
    # them, but alpha, sigma.sq and theta arrive indexed by regime number. left
    # in encounter order they line up only when the first edge starts in regime
    # 1, so half of the sampled histories would be scored with the regimes'
    # parameters permuted.
    state.order <- suppressWarnings(as.numeric(map.states))
    if(!anyNA(state.order)){
      map.states <- map.states[order(state.order)]
    }
  }
  k <- length(map.states)
  tot.states <- factor(map.states)
  tip.states <- factor(data[,1])
  data[,1] <- as.numeric(tip.states)
  #Obtains the state at the root
  root.edge.index <- which(phy$edge[,1] == ntips+1)
  root.state <- which(map.states == names(map[[root.edge.index[2]]][1]))
  ##Begins the construction of the edges matrix -- similar to the ouch format##
  # columns 4 and 5 hold the start and end age of each branch, which weight.mat only
  # reads on the non simmap path or when rescaling. hOUwie is always simmap and never
  # rescales here, so building them would mean a node depth traversal per map per
  # likelihood evaluation for columns nothing goes on to read.
  if(simmap.tree == TRUE & scaleHeight == FALSE){
    if(!is.null(tree.plan)){
      edges <- tree.plan$basic.edges
    }else{
      edges <- cbind(c(1:(n-1)),phy$edge,0,0)
    }
  }else{
    edges <- cbind(c(1:(n-1)),phy$edge,MakeAgeTable(phy, root.age=root.age))
  }
  if(scaleHeight == TRUE){
    Tmax <- max(MakeAgeTable(phy, root.age=root.age))
    edges[,4:5]<-edges[,4:5]/Tmax
    root.age <-  1
    map <- lapply(map, function(x) x/Tmax)
  }
  # column 1 is 1:(n-1) in edge order, so sorting on it undoes the sort on column 3
  # and leaves edges exactly as built. both sorts dropped; edges stays in edge order.

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
  
  Rate.mat <- rbind(alpha, sigma.sq, theta)
  pars <- matrix(c(theta, sigma.sq, alpha), length(theta), 3, dimnames = list(1:length(sigma.sq), c("opt", "sig", "alp")))
  # if the simmap did not simulate every possible state in a given hmm
  if(length(map.states) != dim(Rate.mat)[2]){
    Rate.mat <- Rate.mat[,as.numeric(map.states), drop=FALSE]
    pars <- pars[as.numeric(map.states), , drop=FALSE]
  }
  
  # The dynamic hOUwie path can compute the optimum weights, transformed branch
  # variances, and tip diagonal together while it traverses the map once.
  continuous.moments <- NULL
  if(simmap.tree == TRUE && scaleHeight == FALSE){
    continuous.moments <- continuousMapMoments(
      phy, Rate.mat, pars, root.state = root.state,
      assume.station = !get.root.theta, map = map,
      state.names = map.states,
      edge.order = if(is.null(tree.plan)) NULL else tree.plan$cladewise.order
    )
    W <- continuous.moments$W
  }else if(get.root.theta == TRUE){
    W <- weight.mat(phy, edges, Rate.mat, root.state=root.state,
                    simmap.tree=simmap.tree, root.age=root.age,
                    scaleHeight=scaleHeight, assume.station=FALSE,
                    shift.point=shift.point, map=map,
                    state.names=map.states,
                    edge.order=if(is.null(tree.plan)) NULL else tree.plan$cladewise.order)
  }else{
    W <- weight.mat(phy, edges, Rate.mat, root.state=root.state,
                    simmap.tree=simmap.tree, root.age=root.age,
                    scaleHeight=scaleHeight, assume.station=TRUE,
                    shift.point=shift.point, map=map,
                    state.names=map.states,
                    edge.order=if(is.null(tree.plan)) NULL else tree.plan$cladewise.order)
  }
  
  #Likelihood function for estimating model parameters
  if(get.root.theta == TRUE){
    # maps run root to tip, so only the first name on the root edge is the root state.
    # taking every name would make theta0 as long as that edge has segments and silently
    # recycle the weight matrix against a too long vector of optima.
    root_index <- as.numeric(names(map[[which.min(phy$edge[,1])[1]]])[1])
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
  if(is.null(continuous.moments)){
    transformed.tree <- transformPhy(
      phy, map, pars, tip.paths,
      edge.order=if(is.null(tree.plan)) NULL else tree.plan$cladewise.order
    )
  }else{
    transformed.tree <- continuous.moments[c("tree", "diag")]
  }
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
  nStates <- getNumberOfStates(hOUwie.dat$data.cor[,2])
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
    data.cor <- corProcessData(data.cor, collapse = collapse)
    data.ou <- data.frame(sp = data[,1], 
                          reg = data.cor$corData[,2], 
                          x = data[, dim(data)[2]-1],
                          err = data[, dim(data)[2]])
  }
  if(tip.fog=="none"){
    data.cor <- data[, 1:(dim(data)[2]-1)]
    data.cor <- corProcessData(data.cor)
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

# Sigma squared controls the scale of continuous-trait variation, so its search
# must start on a scale learned from the trait. Using log(2) / tree height here,
# as alpha does, makes the initial stationary variance arbitrary and can send the
# optimizer down the near-Brownian alpha-theta ridge. A tiny positive value keeps
# log-scale optimization defined when the observed trait is constant or contains
# fewer than two finite values; checkStartingUBLB will then raise it to any
# user-supplied lower bound.
getIP.sigma <- function(x){
  finite <- x[is.finite(x)]
  trait.variance <- if(length(finite) > 1) var(finite) else NA_real_
  if(!is.finite(trait.variance) || trait.variance <= 0){
    return(.Machine$double.eps)
  }
  trait.variance
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
# how far a trait has to move right before it can be optimized. the search is over log(p),
# so every parameter it touches - the OU optima included - must be strictly positive, and
# the bounds handed to it are min(trait)/10 and max(trait)*10. landing the shifted minimum
# exactly on zero is therefore no better than leaving it negative, since log(0) is -Inf; it
# has to clear zero by a margin. the margin is a tenth of the trait's own spread rather than
# a tenth of its magnitude because theta is resolved on the log scale: a magnitude-based
# offset would push a tightly clustered trait far from the origin, where neighbouring values
# differ by too little of their log to separate.
getTraitShift <- function(x){
  min_x <- min(x)
  if(min_x >= 0){
    return(0)
  }
  spread <- max(x) - min_x
  # an invariant trait offers no spread to scale the margin by, so use its magnitude, which
  # is nonzero whenever the minimum is negative
  if(spread <= 0){
    spread <- abs(min_x)
  }
  return((0.1 * spread) - min_x)
}
# recover the shift a fit was made under. fits made before the shift was derived from the
# data record only the logical flag, and 50 is the constant they were fit with, so it is the
# only value that puts them back on their original scale
getObjTraitShift <- function(houwie_obj){
  if(!isTRUE(houwie_obj$negative_values)){
    return(0)
  }
  if(is.null(houwie_obj$trait_shift)){
    return(50)
  }
  return(houwie_obj$trait_shift)
}
# organize the houwie output
getHouwieObj <- function(liks_houwie, pars, phy, data, hOUwie.dat, rate.cat, tip.fog, index.disc, index.cont, root.p, nSim, sample_nodes, nStates, discrete_model, continuous_model, time_slice, root.station, get.root.theta,lb_discrete_model,ub_discrete_model,lb_continuous_model,ub_continuous_model,ip, opts, quiet, negative_values, trait_shift=0){
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
    # the thetas occupy the tail of pars, so the count of them has to match the count the
    # caller built pars with - free parameters only, ignoring the cells a model leaves NA
    n_theta <- length(unique(na.omit(index.cont[3,])))
    if(n_theta > 0){
      pars[(length(pars) - n_theta + 1):length(pars)] <- pars[(length(pars) - n_theta + 1):length(pars)] - trait_shift
    }
    if(tip.fog == "none"){
      data[,dim(data)[2]] <- data[,dim(data)[2]] - trait_shift
      solution$solution.ou[3,] <- solution$solution.ou[3,] - trait_shift
    }else{
      data[,dim(data)[2]-1] <- data[,dim(data)[2]-1] - trait_shift
      solution$solution.ou[3,] <- solution$solution.ou[3,] - trait_shift
    }
    # the trait column of data.ou has to come back with the thetas. hOUwie.recon and
    # hOUwie.thorough both pair the returned p against this table, and a theta shifted
    # back against a trait that was not describes a different model entirely.
    hOUwie.dat$data.ou[,3] <- hOUwie.dat$data.ou[,3] - trait_shift
  }
  obj <- list(
    loglik = liks_houwie$TotalLik,
    # effective sample size of the importance weights at the reported optimum.
    # a value near 1 means one sampled history carried the whole likelihood and
    # nSim should be raised (or the proposal improved) before the fit is trusted
    ess = liks_houwie$ess,
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
    sample_nodes = sample_nodes,
    lb_discrete_model=lb_discrete_model, 
    ub_discrete_model=ub_discrete_model,
    lb_continuous_model=lb_continuous_model, 
    ub_continuous_model=ub_continuous_model,
    p=pars,
    ip=ip,
    nSim=nSim,
    opts=opts,
    quiet=quiet,
    negative_values=negative_values,
    trait_shift=trait_shift
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
  # a branch with no state change on it has one segment and nothing to merge. without
  # this 2:length(edge) counts backwards and compares against a name past the end.
  if(length(edge) < 2){
    return(edge)
  }
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


## ---------------------------------------------------------------------------
## Vendored from corHMM
##
## The two functions below are verbatim copies of corHMM's internal
## convertSubHistoryToEdge() (R/makeSimmap.R) and corProcessData()
## (R/corHMM.R), taken from corHMM 2.10.2 (definitions verified unchanged
## through corHMM 2.10.5).
##
## They are copied rather than reached via corHMM::: because CRAN policy
## disallows ::: access to another package's unexported objects. This is a
## TEMPORARY measure: once corHMM exports them, delete these copies and call
## corHMM::convertSubHistoryToEdge() / corHMM::corProcessData() instead.
##
## NOTE: corProcessData() encodes corHMM's data-format contract (the "&"
## polymorphism syntax, "?" for unknown, and the ordering of PossibleTraits).
## If corHMM changes that representation, these copies must be updated in
## lockstep or hOUwie's discrete likelihood will be silently wrong.
##
## corProcessData() calls getRateCatMat(), which corHMM exports and OUwie
## already imports (see NAMESPACE), so no further vendoring is required.
## ---------------------------------------------------------------------------

# convert a substitution history into a mapped edge
convertSubHistoryToEdge <- function(phy, map, RowNames=NULL){
  Traits <- as.numeric(unique(names(unlist(map))))
  if(is.null(RowNames)){
    RowNames <- paste(phy$edge[,1], phy$edge[,2], sep = ",")
  }
  # accumulate each edge's segments into its trait column directly. building the row
  # by row sums over every trait instead costs a closure per trait per edge, and this
  # runs once per simmap for every likelihood evaluation.
  obj <- matrix(0, length(map), length(Traits))
  for(i in seq_along(map)){
    x <- map[[i]]
    col_i <- match(as.numeric(names(x)), Traits)
    for(j in seq_along(x)){
      obj[i, col_i[j]] <- obj[i, col_i[j]] + x[j]
    }
  }
  rownames(obj) <- RowNames
  colnames(obj) <- Traits
  return(obj)
}

corProcessData <- function(data, rate.mat=NULL, collapse=FALSE){
  nCol <- dim(data)[2]
  LevelList <- StateMats <- vector("list", nCol-1)
  # detect the number of states in each column. & is treated as indicating polymorphism. ? is treated as unknown data.
  for(i in 2:nCol){
    data_tmp <- data[,i]
    if(!is.factor(data_tmp)){
      data_tmp <- as.factor(data_tmp)
    }
    States_i <- levels(data_tmp)
    if(any(States_i == "?")){
      States_i <- States_i[!States_i == "?"]
    }
    if(length(grep("&", States_i)) > 0){
      States_i <- unique(unlist(strsplit(States_i, "&")))
    }
    StateMats[[i-1]] <- getRateCatMat(length(States_i))
    LevelList[[i-1]] <- States_i
  }
  # identify the possible trait combinations
  TraitList <- expand.grid(LevelList)
  Traits <- apply(TraitList, 1, function(x) paste(c(x), collapse = "_"))
  # convert each column into a numeric value associated with a member of the trait combinations. ? are associated with all values of that column, & indicates the combination of two or more
  search.strings <- observed.traits_index <- combined.data <- c()
  for(i in 1:dim(data)[1]){
    data_rowi <- data[i,2:nCol]
    # and symbolizes it can be any of the separated states
    search.string_i <- paste("^",paste(sapply(data_rowi, function(x) paste("(", gsub("&", "|", x), ")", sep = "")),collapse = "_"), "$", sep="")
    # ? means it can be any of the states in that character
    search.string_i <- gsub("(?)", ".*", search.string_i, fixed=TRUE)
    # if the data is polymorphic it will now have ands separating the corHMM states
    combined.data[i] <- paste(grep(search.string_i, Traits), collapse="&")
    observed.traits_index <- c(observed.traits_index, grep(search.string_i, Traits))
    search.strings[i] <- search.string_i
  }
  ObservedTraits <- Traits[sort(unique(observed.traits_index))]
  if(collapse){
    corData <- data.frame(sp = data[,1], 
                          d = sapply(search.strings, function(x) 
                            paste(grep(x, ObservedTraits), collapse="&")))
  }else{
    corData <- data.frame(sp = data[, 1], 
                          d = sapply(search.strings, function(x) 
                            paste(grep(x, Traits),collapse = "&")))
  }
  return(list(StateMats = StateMats,  PossibleTraits = Traits, ObservedTraits = ObservedTraits, corData = corData))
}
