
### Ho and Ane test of identifiability of selection optima

#written by Jeremy M. Beaulieu

CheckIdentify <- function(phy, data, simmap.tree=FALSE, get.penalty=TRUE, verbose=TRUE){
    phy = reorder(phy, "pruningwise")
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- n + 1L
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    el <- phy$edge.length
    data.new <- data.frame(data[,2], data[,3], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]
    regime_labs <- c(data.new[,1], phy$node.label)
    #creates n x n punch card of component connections among tips:
    v <- matrix(0, n + phy$Nnode, n)
    #keeps tracks of the regime shifts in the tree:
    regime_shifts <- list()
    shift.number <- 0
    for (i in 1:N) {
        if(des[i] <= n){
            v[des[i],des[i]] = 1
        }
        v[anc[i],] = v[anc[i],]+v[des[i],]
        if(simmap.tree==TRUE){
            print("oops. Forgot!")
        }else{
            if(regime_labs[anc[i]] != regime_labs[des[i]]){
                shift.number <- shift.number + 1
                regime_shifts[[shift.number]] <- i
            }
        }
    }
    identifiable <- check(regime_shifts)
    if(identifiable[1] == 0){
        if(verbose=TRUE){
            cat("The regimes are unidentifiable.", "\n")
        }
        return(identifiable[1])
    }else{
        if(get.penalty == TRUE){
            if(verbose=TRUE){
                cat("The regimes are identifiable.", "\n")
            }
            return(identifiable[2])
        }else{
            if(verbose=TRUE){
                cat("The regimes are identifiable.", "\n")
            }
            return(identifiable[1])
        }
    }
}


##Check function taken from phylolm, which is an embedded function
##for checking identifiability of the OUshifts model:
check <- function(model) {
    ### return identifiability check and modified BIC penalty
    if (length(model)==0) return(c(1,log(n)))
    checkpar = rep(ROOT,n)
    for (i in N:1)
    if (i %in% model) checkpar[which(v[des[i],]==1)] = i
    if (length(levels(as.factor(checkpar))) == length(model)+1) {
        numchange = length(model)
        pen = 3*log(n) + (2*numchange - 1)*log(n)
        for (j in 1:numchange)
        pen = pen + log(length(which(as.factor(checkpar)==levels(as.factor(checkpar))[j])))
        return(c(1,pen))
    } else return(c(0,0))
}

