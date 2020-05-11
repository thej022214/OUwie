#Weight matrix generator taken from Butler and King (2004) and modified to allow multiple alpha parameters

#written by Jeremy M. Beaulieu

weight.mat<-function(phy, edges, Rate.mat, root.state, simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE, assume.station=TRUE, shift.point=0.5){
    
    n <- max(phy$edge[,1])
    ntips <- length(phy$tip.label)
    if(is.null(root.state)) {
        root.state<-which(edges[dim(edges)[1],]==1)-5
        edges <- edges[-1*dim(edges)[1],]
    }
    if(simmap.tree==TRUE){
        k <- length(colnames(phy$mapped.edge))
    }
    if(simmap.tree==FALSE){
        mm <- dim(edges)
        k <- length(6:mm[2])
    }
    pp <- prop.part(phy)
    edges <- edges
    nodevar.root.tot <- rep(0,max(edges[,3]))
    nodevar.k <- rep(0,max(edges[,3]))
    alpha <- Rate.mat[1,]
    W <- matrix(0, ntips, k)

    for(j in 1:k){
        oldregime <- root.state
        n.cov.root.tot = matrix(0, n, 1)
        n.cov.k <- matrix(0, n, 1)
        #Weights for each species per regime
        for(i in 1:length(edges[,1])){
            anc <- edges[i, 2]
            oldtime <- edges[i,4]
            newtime <- edges[i,5]
            if(simmap.tree == TRUE){
                if(scaleHeight == TRUE){
                    currentmap <- phy$maps[[i]]/max(MakeAgeTable(phy, root.age=root.age))
                }
                else{
                    currentmap <- phy$maps[[i]]
                }
            }
            if(simmap.tree==TRUE){
                nodevar.k[i] <- 0
                nodevar.root.tot[i] <- 0
                for (regimeindex in 1:length(currentmap)){
                    regimeduration <- currentmap[regimeindex]
                    newtime <- oldtime + regimeduration
                    regimenumber <- which(colnames(phy$mapped.edge)==names(currentmap)[regimeindex])
                    if(regimenumber == j){
                        nodevar.root.tot[i] <- -alpha[regimenumber]*(newtime-oldtime)
                        nodevar.k[i] <- exp(alpha[regimenumber]*newtime)-exp(alpha[regimenumber]*oldtime)
                    }else{
                        nodevar.root.tot[i] <- -alpha[regimenumber]*(newtime-oldtime)
                        nodevar.k[i] <- nodevar.k[i] + 0
                    }
                    oldtime <- newtime
                }
            }
            if(simmap.tree==FALSE){
                if(anc%in%edges[,3]){
                    start <- which(edges[,3]==anc)
                    oldregime <- which(edges[start,6:(k+5)]==1)
                }
                else{
                    #For the root:
                    oldregime <- root.state
                }
                newregime <- which(edges[i,6:(k+5)]==1)
                if(oldregime==newregime){
                    if(oldregime == j){
                        nodevar.root.tot[i] <- -alpha[oldregime]*(newtime-oldtime)
                        nodevar.k[i] <- exp(alpha[oldregime]*newtime)-exp(alpha[oldregime]*oldtime)
                    }
                    else{
                        nodevar.root.tot[i] <- -alpha[oldregime]*(newtime-oldtime)
                        nodevar.k[i] <- 0
                    }
                }
                else{
                    shifttime <- newtime-((newtime-oldtime) * shift.point)
                    epoch1a <- -alpha[oldregime]*(shifttime-oldtime)
                    epoch1b <- exp(alpha[oldregime]*shifttime)-exp(alpha[oldregime]*oldtime)
                    oldtime <- shifttime
                    epoch2a <- -alpha[newregime]*(newtime-oldtime)
                    epoch2b <- exp(alpha[newregime]*newtime)-exp(alpha[newregime]*oldtime)
                    nodevar.root.tot[i] <- epoch1a + epoch2a
                    
                    if(oldregime==j){
                        nodevar.k[i] <- epoch1b
                    }
                    if(newregime==j){
                        nodevar.k[i] <- epoch2b
                    }
                    if(!newregime==j && !oldregime==j){
                        nodevar.k[i] <- 0
                    }
                }
            }
            n.cov.k[edges[i,3],] <- nodevar.k[i]
            n.cov.root.tot[edges[i,3],] <- nodevar.root.tot[i]
        }
        w.k <- mat.gen(phy, n.cov.k, pp)
        w.root.tot <- mat.gen(phy, n.cov.root.tot, pp)

        W[1:(ntips),j] <- exp(diag(w.root.tot)) * diag(w.k)
    }

    if(assume.station == TRUE){
        w.root.tot <- mat.gen(phy, n.cov.root.tot, pp)
        W[,root.state] <- W[,root.state] + exp(diag(w.root.tot))
    }else{
        w.root.tot <- mat.gen(phy, n.cov.root.tot, pp)
        W <- cbind(exp(diag(w.root.tot)), W)
    }

    #Restandardizes W so that the rows sum to 1 -- Generalized. Will reduce to the simpler model if assuming 1 alpha parameter, but when alpha varies by regime they will sum to 1 (though proportionally should be ok).
    W <- W/rowSums(W)
    W
}


##Matrix generating function taken from vcv.phylo in ape:
mat.gen<-function(phy,piece.wise,pp){
    phy <- reorder(phy, "pruningwise")
    n <- length(phy$tip.label)
    anc <- phy$edge[,1]
    des <- phy$edge[,2]
    ep <- piece.wise[,1]
    comp <- numeric(n + phy$Nnode)
    mat <- matrix(0, n, n)
    
    for (i in length(anc):1) {
        focal <- comp[anc[i]]
        comp[des[i]] <- focal + ep[des[i]]
        j <- i - 1L
        while (anc[j] == anc[i] && j > 0) {
            left <- if (des[j] > n) pp[[des[j] - n]] else des[j]
            right <- if (des[i] > n) pp[[des[i] - n]] else des[i]
            mat[left, right] <- mat[right, left] <- focal
            j <- j - 1L
        }
    }
    diag.elts <- 1 + 0:(n - 1)*(n + 1)
    mat[diag.elts] <- comp[1:n]
    
    mat
}

