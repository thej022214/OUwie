#OU functions for using three-point algorithm of Ho and Ane 2013

#written by James D. Boyko


# takes a node based reconstruction and transforms it into a simmap
makeSimmapFromNode <- function(phy){
    Map <- vector("list", dim(phy$edge)[1])
    NodeLabels <- phy$node.label[phy$edge[,1] - Ntip(phy)]
    for(i in 1:dim(phy$edge)[1]){
        tmp <- phy$edge.length[i]
        names(tmp) <- NodeLabels[i]
        Map[[i]] <- tmp
    }
    phy$maps <- Map
    phy$mapped.edge <- corHMM:::convertSubHistoryToEdge(tree, Map)
    class(phy) <-  c("simmap", setdiff(class(phy), "simmap"))
    return(phy)
}


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


# transforms the phylogeny based on a set of paramaters and a simmap
transformPhy <- function(phy, pars){
    # phy must be of class simmap
    nTip <- length(phy$tip.label)
    RootAge <- max(branching.times(phy))
    NodeAges <- branching.times(phy)[phy$edge[,1] - nTip]
    ModMap <- Map <- phy$maps
    D <- V_Tilde <- numeric(dim(phy$edge)[1])
    for(i in 1:dim(phy$edge)[1]){
        # evaluate the map for this particular edge and calculate the tipward variance
        NodeAge_i <- NodeAges[i]
        DistRoot_i <- RootAge - NodeAge_i
        Map_i <- Map[[i]]
        # the age of epoch j starts at the node age
        # EpochAge_j <- NodeAge_i
        Dist_rootward <- DistRoot_i
        w <- v <- 0
        for(j in 1:length(Map_i)){
            # distance the root of epoch j starts at the node distance and ends at node dist + epoch length
            Dist_tipward <- Dist_rootward + Map_i[j]
            # the length of the epoch is scaled by the alpha parameter of that epoch
            Sigma_j <- pars[,2][match(names(Map_i)[j], rownames(pars))]
            Alpha_j <- pars[,3][match(names(Map_i)[j], rownames(pars))]
            # calculate the descendent distance from the root based on a fixed root distribution
            tmp <- Sigma_j * (exp(2 * Alpha_j * Dist_tipward) - exp(2 * Alpha_j * Dist_rootward))/2/Alpha_j
            v <- v + tmp
            ModMap[[i]][j] <- tmp
            w <- w + (Alpha_j * (Dist_tipward - Dist_rootward))
            # The new distance from nodes
            Dist_rootward <- Dist_tipward
        }
        V_Tilde[i] <- v
        D[i] <- w
    }
    
    # calculates the diagonal matrix for each tip i
    DiagWt <- numeric(nTip)
    names(DiagWt) <- phy$tip.label
    for(i in 1:nTip){
        DiagWt[i] <- exp(-sum(D[getPathToRoot(phy, i)]))
    }
    
    phy$edge.length <- V_Tilde
    phy$maps <- ModMap
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





