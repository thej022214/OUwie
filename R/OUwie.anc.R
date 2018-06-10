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


anc.likelihood <- function(x, fitted.OUwie.object) {
    fitted.OUwie.object$data[grepl("node_", fitted.OUwie.object$data[,1]),3] <- x
    return(OUwie.fixed(fitted.OUwie.object$phy,fitted.OUwie.object$data, fitted.OUwie.object$model,fitted.OUwie.object$simmap.tree, fitted.OUwie.object$root.age, fitted.OUwie.object$scaleHeight,fitted.OUwie.object$root.station, fitted.OUwie.object$alpha, fitted.OUwie.object$sigma.sq, fitted.OUwie.object$theta, fitted.OUwie.object$clade, fitted.OUwie.object$mserr, quiet=TRUE))
}

OUwie.anc<-function(fitted.OUwie.object){
    nloptr(x0=rep(median(fitted.OUwie.object$data[,3], na.rm=TRUE), sum(grepl("node_", data[,1]))), eval_f=likelihood, _____)
}
