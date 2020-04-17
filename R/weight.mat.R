#Weight matrix generator taken from Butler and King (2004) and modified to allow multiple alpha parameters

#written by Jeremy M. Beaulieu

weight.mat<-function(phy, edges, Rate.mat, root.state, simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE, assume.station=TRUE, shift.point=.5){
    
    n=max(phy$edge[,1])
    ntips=length(phy$tip.label)
    if(is.null(root.state)) {
        root.state<-which(edges[dim(edges)[1],]==1)-5
        edges<-edges[-1*dim(edges)[1],]
    }
    if(simmap.tree==TRUE){
        k=length(colnames(phy$mapped.edge))
    }
    if(simmap.tree==FALSE){
        mm<-dim(edges)
        k<-length(6:mm[2])
    }
    pp <- prop.part(phy)
    nodevar=rep(0,max(edges[,3]))
    alpha=Rate.mat[1,]
    if(assume.station==TRUE){
        W <- matrix(0,ntips,k)
        oldregime=root.state
        for(j in 1:k){
            n.cov=matrix(0, n, 1)
            #Weights for each species per regime
            for(i in 1:length(edges[,1])){
                anc = edges[i, 2]
                oldtime = edges[i,4]
                newtime = edges[i,5]
                if(simmap.tree==TRUE){
                    if(scaleHeight==TRUE){
                        currentmap <- phy$maps[[i]]/max(MakeAgeTable(phy, root.age=root.age))
                    }
                    else{
                        currentmap <- phy$maps[[i]]
                    }
                }
                if(simmap.tree==TRUE){
                    nodevar[i] <- 0
                    for (regimeindex in 1:length(currentmap)){
                        regimeduration <- currentmap[regimeindex]
                        newtime <- oldtime + regimeduration
                        regimenumber <- which(colnames(phy$mapped.edge)==names(currentmap)[regimeindex])
                        if(regimenumber == j){
                            nodevar[i] <- exp(-alpha[root.state])*(exp(alpha[regimenumber]*newtime)-exp(alpha[regimenumber]*oldtime))
                        }
                        else{
                            nodevar[i] <- nodevar[i]
                        }
                        oldtime <- newtime
                    }
                }
                if(simmap.tree==FALSE){
                    if(anc%in%edges[,3]){
                        start=which(edges[,3]==anc)
                        oldregime=which(edges[start,6:(k+5)]==1)
                    }
                    else{
                        #For the root:
                        oldregime=root.state
                    }
                    newregime=which(edges[i,6:(k+5)]==1)
                    if(oldregime==newregime){
                        if(newregime == j){
                            nodevar[i] <- exp(-alpha[root.state])*(exp(alpha[oldregime]*newtime)-exp(alpha[oldregime]*oldtime))
                        }
                        else{
                            nodevar[i] <- 0
                        }
                    }
                    else{
                        shifttime=newtime-((newtime-oldtime)*shift.point)
                        epoch1 = exp(-alpha[root.state])*(exp(alpha[oldregime]*shifttime)-exp(alpha[oldregime]*oldtime))
                        oldtime = shifttime
                        newtime = newtime
                        epoch2 = exp(-alpha[root.state])*(exp(alpha[newregime]*newtime)-exp(alpha[newregime]*oldtime))
                        if(oldregime==j){
                            nodevar[i]=epoch1
                        }
                        if(newregime==j){
                            nodevar[i]=epoch2
                        }
                        if(!newregime==j && !oldregime==j){
                            nodevar[i] = 0
                        }
                    }
                }
                n.cov[edges[i,3],] = nodevar[i]
            }
            w.piece <- mat.gen(phy, n.cov, pp)
            W[1:(ntips),j] <- diag(w.piece)
        }
        W[,root.state] <- W[,root.state]+exp(-alpha[root.state])
        #W <- W/rowSums(W)
    }
    
    if(assume.station==FALSE){
        W <- matrix(0,ntips,k+1)
        W.piece.root <- matrix(0, ntips, ntips)
        for(j in 1:k){
            oldregime=root.state
            n.cov=matrix(0, n, 1)
            #Weights for each species per regime
            for(i in 1:length(edges[,1])){
                anc = edges[i, 2]
                oldtime=edges[i,4]
                newtime=edges[i,5]
                if(simmap.tree==TRUE){
                    if(scaleHeight==TRUE){
                        currentmap<-phy$maps[[i]]/max(MakeAgeTable(phy, root.age=root.age))
                    }
                    else{
                        currentmap<-phy$maps[[i]]
                    }
                }
                if(simmap.tree==TRUE){
                    nodevar1[i]=0
                    nodevar2[i]=0
                    for (regimeindex in 1:length(currentmap)){
                        regimeduration <- currentmap[regimeindex]
                        newtime <- oldtime + regimeduration
                        regimenumber <- which(colnames(phy$mapped.edge)==names(currentmap)[regimeindex])
                        if (regimenumber == j){
                            nodevar[i] <- exp(-alpha[root.state])*(exp(alpha[regimenumber]*newtime)-exp(alpha[regimenumber]*oldtime))
                        }
                        else{
                            nodevar[i] <- nodevar1[i]
                        }
                        oldtime <- newtime
                    }
                }
                if(simmap.tree==FALSE){
                    if(anc%in%edges[,3]){
                        start=which(edges[,3]==anc)
                        oldregime=which(edges[start,6:(k+5)]==1)
                    }
                    else{
                        #For the root:
                        oldregime=root.state
                    }
                    newregime=which(edges[i,6:(k+5)]==1)
                    if(oldregime==newregime){
                        if(newregime==j){
                            nodevar[i]=exp(-alpha[root.state])*(exp(alpha[regimenumber]*newtime)-exp(alpha[regimenumber]*oldtime))
                        }
                        else{
                            nodevar[i] <- 0
                        }
                    }
                    else{
                        shifttime=newtime-((newtime-oldtime)*shift.point)
                        epoch1 = exp(-alpha[root.state])*(exp(alpha[oldregime]*shifttime)-exp(alpha[oldregime]*oldtime))
                        oldtime = shifttime
                        newtime = newtime
                        epoch2 = exp(-alpha[root.state])*(exp(alpha[newregime]*newtime)-exp(alpha[newregime]*oldtime))
                        if(oldregime==j){
                            nodevar[i]=epoch1
                        }
                        if(newregime==j){
                            nodevar[i]=epoch2
                        }
                        if(!newregime==j & !oldregime==j){
                            nodevar[i]=0
                        }
                    }
                }
                n.cov[edges[i,3]]=nodevar[i]
            }
            w.piece <- mat.gen(phy,n.cov,pp)
            W[1:(ntips),j+1] <- diag(w.piece)
        }
        W[,root.state] <- W[,root.state]+exp(-alpha[root.state])
        #W <- W/rowSums(W)
    }
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
