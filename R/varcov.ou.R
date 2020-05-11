#OU variance-covariance matrix generator

#written by Jeremy M. Beaulieu

varcov.ou <- function(phy, edges, Rate.mat, root.state, simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE, assume.station=TRUE, shift.point=.5){
    
    if(assume.station == TRUE){
        alpha=Rate.mat[1,1]
        sigma=Rate.mat[2,1]
        vcv <- quickVCV(phy=phy, alpha=alpha, sigma.sq=sigma, scaleHeight=scaleHeight)
    }else{
        if(is.null(root.state)) {
            root.state<-which(edges[dim(edges)[1],]==1)-5
            edges<-edges[-1*dim(edges)[1],]
        }
        n=max(phy$edge[,1])
        ntips=length(phy$tip.label)
        if(simmap.tree==TRUE){
            k=length(colnames(phy$mapped.edge))
        }
        if(simmap.tree==FALSE){
            mm<-dim(edges)
            k<-length(6:mm[2])
        }
        pp <- prop.part(phy)
        oldregime=root.state
        nodevar1=rep(0,max(edges[,3]))
        nodevar2=rep(0,max(edges[,3]))
        alpha=Rate.mat[1,]
        sigma=Rate.mat[2,]
        n.cov1=matrix(rep(0,n), n, 1)
        n.cov2=matrix(rep(0,n), n, 1)
        
        if(simmap.tree==TRUE){
            regimeindex<-colnames(phy$mapped.edge)
            for(i in 1:length(edges[,1])){
                anc = edges[i, 2]
                desc = edges[i, 3]
                
                if(scaleHeight==TRUE){
                    currentmap<-phy$maps[[i]]/max(MakeAgeTable(phy, root.age=root.age))
                }
                else{
                    currentmap <- phy$maps[[i]]
                }
                oldtime=edges[i,4]
                for (regimeindex in 1:length(currentmap)){
                    regimeduration <- currentmap[regimeindex]
                    newtime <- oldtime+regimeduration
                    regimenumber <- which(colnames(phy$mapped.edge)==names(currentmap)[regimeindex])
                    nodevar1[i] <- nodevar1[i]+alpha[regimenumber]*(newtime-oldtime)
                    nodevar2[i] <- nodevar2[i]+sigma[regimenumber]*((exp(2*alpha[regimenumber]*newtime)-exp(2*alpha[regimenumber]*oldtime))/(2*alpha[regimenumber]))
                    oldtime <- newtime
                    newregime <- regimenumber
                }
                oldregime=newregime
                n.cov1[edges[i,3],]=nodevar1[i]
                n.cov2[edges[i,3],]=nodevar2[i]
            }
        }
        if(simmap.tree==FALSE){
            for(i in 1:length(edges[,1])){
                anc <- edges[i,2]
                oldtime <- edges[i,4]
                newtime <- edges[i,5]
                if(anc%in%edges[,3]){
                    start <- which(edges[,3]==anc)
                    oldregime <- which(edges[start,6:(k+5)]==1)
                }
                else{
                    #For the root:
                    oldregime=root.state
                }
                newregime=which(edges[i,6:(k+5)]==1)
                if(oldregime==newregime){
                    nodevar1[i] <- alpha[oldregime]*(newtime-oldtime)
                    nodevar2[i] <- sigma[oldregime]*((exp(2*alpha[oldregime]*newtime)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))
                }
                else{
                    shifttime <- newtime-((newtime-oldtime)*shift.point)
                    epoch1a <- alpha[oldregime]*(shifttime-oldtime)
                    epoch1b <- sigma[oldregime]*((exp(2*alpha[oldregime]*shifttime)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))
                    oldtime <- shifttime
                    newtime <- newtime
                    epoch2a <- alpha[newregime]*(newtime-oldtime)
                    epoch2b <- sigma[newregime]*((exp(2*alpha[newregime]*newtime)-exp(2*alpha[newregime]*oldtime))/(2*alpha[newregime]))
                    nodevar1[i] <- epoch1a+epoch2a
                    nodevar2[i] <- epoch1b+epoch2b
                }
                oldregime <- newregime
                n.cov1[edges[i,3],] <- nodevar1[i]
                n.cov2[edges[i,3],] <- nodevar2[i]
            }
        }
        vcv1 <- mat.gen(phy,n.cov1,pp)
        vcv2 <- mat.gen(phy,n.cov2,pp)
        if(any(abs(diff(alpha)) > 0)){
            species.variances <- diag(vcv1)
            species.total.variances <- matrix(0, dim(vcv1)[2], dim(vcv1)[2])
            count=0
            for(i in 1:dim(vcv1)[2]) {
                for(j in 1:dim(vcv1)[2]){
                    if(i==j){
                        break;
                    }else{
                        species.total.variances[i,j] <- exp(-(species.variances[i] + species.variances[j]))
                        count=count+1
                    }
                }
            }
            diag(species.total.variances) <- exp(-2 * diag(vcv1))
            vcv <- species.total.variances * vcv2
        }else{
            if(is.null(root.age)){
                root.age <- max(branching.times(phy))
            }
            vcv <- exp(-2*alpha[1]*max(root.age)) * vcv2
        }
    }
    vcv
}


## Quick VCV maker of OU1 and OUM -- since alpha and sigma.sq are constants, and since the regime does not matter, we can just do a simple plug and chug. It's a tad slower, but mostly for testing purposes.
quickVCV <- function(phy, alpha, sigma.sq, scaleHeight){
    phy$node.label <- NULL
    vcv <- matrix(0, Ntip(phy), Ntip(phy))
    split.times <- branching.times(phy)
    if(scaleHeight == TRUE){
        split.times <- split.times/max(split.times)
    }
    tot.time <- max(split.times)
    for(i in 1:Ntip(phy)){
        for(j in 1:Ntip(phy)){
            if(i == j){
                dij <- 0
                tij <- tot.time
                vcv[i,j] <- (sigma.sq /(2 * alpha)) * exp(-2 * alpha * dij)
            }else{
                split <- getMRCA(phy, tip=c(phy$tip.label[i], phy$tip.label[j]))
                dij <- split.times[which(names(split.times)==split)]
                tij <- tot.time - dij
                vcv[i,j] <- (sigma.sq /(2 * alpha)) * exp(-2 * alpha * dij)
            }
        }
    }
    return(vcv)
}


##Matrix generating function taken from vcv.phylo in ape:
mat.gen <- function(phy,piece.wise,pp){
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

