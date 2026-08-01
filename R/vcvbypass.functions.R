#OU functions for using three-point algorithm of Ho and Ane 2013

#written by James D. Boyko


## takes a node based reconstruction and returns a map (identical to a map from simmap)
getMapFromNode <- function(phy, tipstates, nodestates, shift.point){
  Map <- vector("list", dim(phy$edge)[1])
  Data <- c(tipstates, nodestates)
  NodeStates <- cbind(Data[phy$edge[,1]], Data[phy$edge[,2]])
  for(i in 1:dim(phy$edge)[1]){
    from <- as.character(NodeStates[i,1])
    to <- as.character(NodeStates[i,2])
    if(from == to){
      tmp <- phy$edge.length[i]
      names(tmp) <- from
      Map[[i]] <- tmp
    }else{
      shift.time <- shift.point * phy$edge.length[i]
      tmp <- c(phy$edge.length[i] - shift.time, shift.time)
      names(tmp) <- c(from, to)
      Map[[i]] <- tmp
    }
  }
  return(Map)
}

# data(tworegime)
# phy <- tree
# tipstates <- trait[,2]
# nodestates <- phy$node.label
# shift.point <- 0.5
# getMapFromNode(phy, round(runif(length(phy$tip.label), 1, 2)), nodestates, 0.5)
# getMapFromNode(phy, tipstates, nodestates, 0)


# gets the path from a vertex to the root as an index of the edge matrix
getPathToRoot <- function(phy, tip){
  nTip <- length(phy$tip.label)
  root <- nTip + 1
  path <- 0
  count <- 1
  while(tip != root){
    tip.ind <- which(phy$edge[,2] == tip)
    path <- c(path, tip.ind)
    count <- count + 1
    tip <- phy$edge[tip.ind,1]
  }
  path <- path[-1]
  
  return(path)
}


# Compute all continuous-map quantities used by OUwie.basic in one cladewise
# traversal. This is the deep internal interface used by dynamic hOUwie: callers
# provide a tree and an edge-segment map and receive the optimum weights, transformed
# branch variances, and tip attenuation together.
continuousMapMoments <- function(phy, Rate.mat, pars, root.state,
		assume.station = TRUE, map = NULL, state.names = NULL,
		edge.order = NULL) {
	nTip <- length(phy$tip.label)
	nState <- ncol(Rate.mat)
	if(is.null(map)) map <- phy$maps
	if(is.null(state.names)) state.names <- colnames(phy$mapped.edge)
	state.index <- stats::setNames(seq_len(nState), state.names)
	weight.alpha <- Rate.mat[1, ]
	transform.alpha <- pars[, 3]
	transform.sigma.sq <- pars[, 2]

	if(is.null(edge.order)){
		edge.order <- ape::reorder.phylo(phy, "cladewise", index.only = TRUE)
	}
	nNodeTotal <- max(phy$edge)

	# Ato is the integral of alpha from the root. cum.weight stores, for all
	# regimes at once, the unattenuated optimum contribution along that path.
	Ato.weight <- numeric(nNodeTotal)
	Ato.transform <- numeric(nNodeTotal)
	cum.weight <- matrix(0, nNodeTotal, nState)
	transformed.variance <- numeric(nrow(phy$edge))

	for(edge.index in edge.order){
		ancestor <- phy$edge[edge.index, 1]
		descendant <- phy$edge[edge.index, 2]
		current.map <- map[[edge.index]]
		regimes <- unname(state.index[names(current.map)])
		transform.regimes <- match(names(current.map), rownames(pars))
		durations <- as.numeric(current.map)
		current.weight.alpha <- Ato.weight[ancestor]
		current.transform.alpha <- Ato.transform[ancestor]
		edge.weight <- numeric(nState)
		variance <- 0

		for(segment.index in seq_along(durations)){
			regime <- regimes[segment.index]
			duration <- durations[segment.index]
			transform.regime <- transform.regimes[segment.index]
			segment.weight.alpha <- weight.alpha[regime]
			segment.transform.alpha <- transform.alpha[transform.regime]
			segment.sigma.sq <- transform.sigma.sq[transform.regime]
			end.weight.alpha <- current.weight.alpha +
				segment.weight.alpha * duration
			end.transform.alpha <- current.transform.alpha +
				segment.transform.alpha * duration

			edge.weight[regime] <- edge.weight[regime] +
				(exp(end.weight.alpha) - exp(current.weight.alpha))
			if(segment.transform.alpha == 0){
				variance <- variance + segment.sigma.sq *
					exp(2 * current.transform.alpha) * duration
			}else{
				variance <- variance + segment.sigma.sq *
					(exp(2 * end.transform.alpha) -
					 exp(2 * current.transform.alpha)) /
					(2 * segment.transform.alpha)
			}
			current.weight.alpha <- end.weight.alpha
			current.transform.alpha <- end.transform.alpha
		}

		Ato.weight[descendant] <- current.weight.alpha
		Ato.transform[descendant] <- current.transform.alpha
		cum.weight[descendant, ] <- cum.weight[ancestor, ] + edge.weight
		transformed.variance[edge.index] <- variance
	}

	root.attenuation <- exp(-Ato.weight[seq_len(nTip)])
	W <- cum.weight[seq_len(nTip), , drop = FALSE] * root.attenuation
	if(assume.station){
		W[, root.state] <- W[, root.state] + root.attenuation
	}else{
		W <- cbind(root.attenuation, W)
	}
	W <- W / rowSums(W)

	transformed.tree <- phy
	transformed.tree$edge.length <- transformed.variance
	transformed.tree$maps <- NULL
	transformed.tree$mapped.edge <- NULL
	attr(transformed.tree, "map.order") <- NULL
	class(transformed.tree) <- "phylo"

	list(
		W = W,
		tree = transformed.tree,
		diag = stats::setNames(exp(-Ato.transform[seq_len(nTip)]), phy$tip.label)
	)
}


# transforms the phylogeny based on a set of parameters and a simmap
transformPhy <- function(phy, map, pars, tip.paths = NULL, edge.order = NULL) {
	#phy must be of class simmap
	nTip <- length(phy$tip.label)

	#Ato[ancestor] must be filled before the edge below it is visited, so edges are
	#walked in cladewise order. map and tip.paths are indexed against the tree as the
	#caller supplied it, so we permute rather than reorder phy and undo it below.
	ord <- edge.order
	if(is.null(ord)){
		ord <- ape::reorder.phylo(phy, "cladewise", index.only = TRUE)
	}
	edge <- phy$edge[ord, , drop = FALSE]

	nEdge <- nrow(edge)
	nNodeTot <- nTip + phy$Nnode

	Ato <- numeric(nNodeTot)
	root <- nTip + 1
	Ato[root] <- 0
	D <- V_Tilde <- numeric(nEdge)

	for (i in 1:nEdge) {
		anc <- edge[i, 1]
		des <- edge[i, 2]
		Map_i <- map[[ord[i]]]
		Acur <- Ato[anc]
		w <- 0
		v <- 0
		for (j in 1:length(Map_i)) {
			dt <- as.numeric(Map_i[j])
			Sigma_j <- pars[, 2][match(names(Map_i)[j], rownames(pars))]
			Alpha_j <- pars[, 3][match(names(Map_i)[j], rownames(pars))]
			tmp.w <- Alpha_j * dt
			w <- w + tmp.w
			if (Alpha_j == 0) {
				tmp.v <- Sigma_j * exp(2 * Acur) * dt
			}else{
				Aend <- Acur + Alpha_j * dt
				tmp.v <- Sigma_j * (exp(2 * Aend) - exp(2 * Acur)) / (2 * Alpha_j)
			}
			v <- v + tmp.v
			Acur <- Acur + Alpha_j * dt
		}
		V_Tilde[ord[i]] <- v
		D[ord[i]] <- w
		Ato[des] <- Acur
	}
	phy$edge.length <- V_Tilde
	#calculates the diagonal matrix for each tip i
	DiagWt <- numeric(nTip)
	names(DiagWt) <- phy$tip.label

	if (is.null(tip.paths)) {
		for (i in 1:nTip) {
			DiagWt[i] <- exp(-sum(D[getPathToRoot(phy, i)]))
		}
	}else{
		for (i in 1:nTip) {
			DiagWt[i] <- exp(-sum(D[tip.paths[[i]]]))
		}
	}
	#the transformed edge lengths no longer agree with the regime painting, so maps and
	#mapped.edge are dropped rather than carried along stale. this also keeps the tree a
	#plain phylo: three.point.compute reorders what it is handed, and on a simmap that
	#dispatches to reorderSimmap, which permutes those two elements for nothing.
	phy$maps <- NULL
	phy$mapped.edge <- NULL
	attr(phy, "map.order") <- NULL
	class(phy) <- "phylo"

	obj <- list(tree = phy, diag = DiagWt)
	return(obj)
}


# transforms the phylogeny based on a set of parameters and a simmap
transformPhy.old <- function(phy, map, pars, tip.paths=NULL){
  # phy must be of class simmap
  nTip <- length(phy$tip.label)
  RootAge <- max(branching.times(phy))
  NodeAges <- branching.times(phy)[phy$edge[,1] - nTip]
  # ModMap <- Map <- map
  D <- V_Tilde <- numeric(dim(phy$edge)[1])
  for(i in 1:dim(phy$edge)[1]){
	# evaluate the map for this particular edge and calculate the tipward variance
	NodeAge_i <- NodeAges[i]
	DistRoot_i <- RootAge - NodeAge_i
	Map_i <- map[[i]]
	# the age of epoch j starts at the node age
	Dist_rootward <- DistRoot_i
	z <- w <- v <- 0
	for(j in 1:length(Map_i)){
	  # distance the root of epoch j starts at the node distance and ends at node dist + epoch length
	  Dist_tipward <- Dist_rootward + Map_i[j]
	  # the length of the epoch is scaled by the alpha parameter of that epoch
	  Sigma_j <- pars[,2][match(names(Map_i)[j], rownames(pars))]
	  Alpha_j <- pars[,3][match(names(Map_i)[j], rownames(pars))]
	  # calculate the descendant distance from the root based on a fixed root distribution
	  tmp.w <- Alpha_j * (Dist_tipward - Dist_rootward)
	  tmp.v <- Sigma_j * (exp(2 * Alpha_j * Dist_tipward) - exp(2 * Alpha_j * Dist_rootward))/2/Alpha_j
	  v <- v + tmp.v
	  w <- w + tmp.w
	  # ModMap[[i]][j] <- tmp.v
	  # The new distance from nodes
	  Dist_rootward <- Dist_tipward
	}
	V_Tilde[i] <- v
	D[i] <- w
  }
  phy$edge.length <- V_Tilde
  # calculates the diagonal matrix for each tip i
  DiagWt <- numeric(nTip)
  names(DiagWt) <- phy$tip.label
  if(is.null(tip.paths)){
	for(i in 1:nTip){
	  DiagWt[i] <- exp(-sum(D[getPathToRoot(phy, i)]))
	}
  }else{
	for(i in 1:nTip){
	  DiagWt[i] <- exp(-sum(D[tip.paths[[i]]]))
	}
  }
  obj <- list(tree = phy, diag = DiagWt)
  return(obj)
}



getOULik <- function(phy, y, X, pars){
  # transform the phylogeny based on params
  tre <- transformPhy(phy, pars)
  # use the transformed phylogeny for the three point algorithm
  comp <- three.point.compute(tre$tree, y, X, tre$diag)
  # calculate the likelihood
  lik <- -as.numeric(Ntip(phy) * log(2 * pi) + comp$logd + (comp$PP - 2 * comp$QP + comp$QQ))/2
  return(lik)
}




