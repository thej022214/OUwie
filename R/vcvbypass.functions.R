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


# transforms the phylogeny based on a set of parameters and a simmap
transformPhy <- function(phy, map, pars, tip.paths=NULL){
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






