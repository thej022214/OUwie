#OUwie likelihood calculator

#written by Jeremy M. Beaulieu and Brian O'Meara

#Allows the user to calculate the likelihood given a specified set of parameter values while estimating the states at internal nodes. Assumes you have estimated the paramters already using OUwie

## NOTE THIS IS NOT WORKING YET.

# Idea is to take a tree, add tips at each node, add tips to the data (including figuring out regimes), and optimize the data for these made up tips. Take in OUwie object

attach.stub.taxon <- function(node, phy, tip.name=NULL) {
  if(is.null(tip.name)) {
    tip.name <- paste0("node_", node)
  }
  return(ape::bind.tree(phy, structure(list(edge = structure(c(2L, 1L), .Dim = 1:2), tip.label = tip.name, Nnode = 1L, edge.length = 0), .Names = c("edge", "tip.label", "Nnode", "edge.length"), class = "phylo"), where=node))
}

attach.stub.taxa <- function(phy) {
    root.node <- ape::Ntip(phy)+1
    start.node <- root.node+1
    n.node.not.root <- ape::Nnode(phy)-1
    for (i in sequence(n.node.not.root)) {
        phy <- attach.stub.taxon(start.node + 2*(i-1), phy, tip.name = paste0("node_", root.node+i)) #all the node numbers go up by one as we add tips
    }
    return(phy)
}

add.stub.taxon.to.data <- function(taxon, phy, data) {
    data[,1] <- as.character(data[,1])
    new.row <- data[1,]
    new.row[1,] <- NA
    new.row[1,1] <- taxon
    if(ncol(new.row)==4) { #we have measurement error
        new.row[1,4] <- 0 #but not in our internal node
    }
    taxon.id <- which(phy$tip.label==taxon)
    focal.node.id <- phy$edge[which(phy$edge[,2]==taxon.id),1]
    new.row[1,2] <- phy$node.label[focal.node.id - ape::Ntip(phy)]
    return(rbind(data, new.row))
}

add.stub.taxa.to.data <- function(phy, data) {
    taxa.to.add <- phy$tip.label[grepl("node_", phy$tip.label)]
    for (i in sequence(length(taxa.to.add))) {
        data <- add.stub.taxon.to.data(taxa.to.add[i], phy, data)
    }
    data <- data[match(phy$tip.label, data[,1]), ]
    return(data)
}

#So, OUwie data going in has a col for taxon names. Coming out of OUwie (and saved in OUwie object), it has rownames instead. So use these functions. Convert the OUwie object to have a regular OUwie formatted input
ouwie.col.to.row.names <- function(x) {
    rownames(x) <- x[,1]
    return(x[,-1])
}

ouwie.row.names.to.col <- function(x) {
    return(cbind(Genus_species=rownames(x), x))
}



anc.likelihood <- function(x, fitted.OUwie.object) {
    traits <- fitted.OUwie.object$data
    traits[grepl("node_", traits[,1]),3] <- x
    #model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"),simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE,root.station=TRUE, alpha=NULL, sigma.sq=NULL, theta=NULL, clade=NULL, mserr="none", quiet=FALSE)
    return(-1*OUwie.fixed(phy=fitted.OUwie.object$phy,data=traits, model=fitted.OUwie.object$model, simmap.tree=fitted.OUwie.object$simmap.tree, root.age=fitted.OUwie.object$root.age, scaleHeight=FALSE, root.station=fitted.OUwie.object$root.station, alpha=fitted.OUwie.object$solution['alpha',], sigma.sq=fitted.OUwie.object$solution['sigma.sq',], theta=fitted.OUwie.object$theta[,1], clade=NULL, mserr=ifelse(is.null(fitted.OUwie.object$mserr.est), "none", fitted.OUwie.object$mserr.est), quiet=TRUE)$loglik)
}

OUwie.anc<-function(fitted.OUwie.object, opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="100000", "ftol_rel"=.Machine$double.eps^0.5)){
    fitted.OUwie.object$phy <- attach.stub.taxa(fitted.OUwie.object$phy)
    fitted.OUwie.object$data <- add.stub.taxa.to.data(fitted.OUwie.object$phy, ouwie.row.names.to.col(fitted.OUwie.object$data))
    traits <- fitted.OUwie.object$data
    result <- nloptr(x0=rep(median(traits[,3], na.rm=TRUE), sum(grepl("node_", traits[,1]))), eval_f=anc.likelihood, opts=opts, fitted.OUwie.object=fitted.OUwie.object)
    traits[grepl("node_", traits[,1]),3] <- result$solution
    fitted.OUwie.object$data <- ouwie.col.to.row.names(traits)
    class(fitted.OUwie.object) <- c("OUwieRecon", "OUwie")
    return(fitted.OUwie.object)
}

plot.OUwieRecon <- function(x, ...) {
    quantitative.trait <- x$data[,2]
    names(quantitative.trait) <- rownames(x$data)
    phytools::contMap(x$phy,quantitative.trait)
}