# Functions to deal with bad matrix condition

# Written by Brian O'Meara

tree.kappa <- function(phy){
    return(kappa(ape::vcv(phy)))
}

fix.kappa <- function(phy, data, threshold = log(40)) {
    current.kappa <- tree.kappa(phy)
    while(current.kappa > threshold & ape::Ntip(phy)>2) {
        tree.structure <- cbind(phy$edge, phy$edge.length)
        tip.structure <- tree.structure[which(tree.structure[,2]<=ape::Ntip(phy)),]
        short.tip <- tip.structure[which.min(tip.structure[,3]),2]
        data <- data[data[,1]!=phy$tip.label[short.tip],]
        phy <- ape::drop.tip(phy, short.tip)
        current.kappa <- tree.kappa(phy)
    }
    return(list(phy=phy, data=data))
}
