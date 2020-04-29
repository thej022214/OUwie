library(phytools)
#need add getDescendants to the NAMESPACE

#criterion: If aicc or aic, use rgenoud, where the fitness is the score (and try to minimize)
#alpha.max.k=3: allows for three OU regimes and one regime where alpha is set to ~0 (brownian motion)
#    if alpha.max.k=0, only BM models are investigated
#sigma.max.k, theta.max.k: work in same way, but must be at least 1 in each case
#note that these params can be stuck together in interesting ways. For example, a BM jump model has multiple thetas but constant everything else
#There are as many free parameters for rgenoud as the number of nodes (not edges) * 3: params for theta, sigma, and alpha
#so each individual for rgenoud has theta1, theta2, theta3....,sigma1, sigma2, sigma3, ....,alpha1, alpha2, alpha3, where theta1 is the mapping to the theta free parameter for node 1, and so forth

OUwie.dredge <- function(phy, data, criterion=c("AIC", "AICc", "BIC", "mBIC"), shift.max=3, sigma.max.k=3, alpha.max.k=3, root.station=FALSE, shift.point=0.5, opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5)) {
    
    #Coerce the data so that it will run in OUwie -- using values of OU1 as the starting points:
    cat("Initializing...","\n")
    
    data2 <- data.frame(taxon=as.character(data[,1]), regime=sample(c(1:2), length(data[,1]), replace=TRUE), trait=data[,2], stringsAsFactors=FALSE)
    phy$node.label <- sample(c(1:2),phy$Nnode, replace=TRUE)
    start.vals <- OUwie(tree, data2, model=c("OU1"), quiet=TRUE, clade=c("t1","t2"), root.station=FALSE)

    cat("Begin optimization routine -- Starting values:", c(start.vals$solution[2,1], start.vals$solution[1,2]), "\n")
    
    find.shifts <- GetShiftModel(phy, data, nmax=shift.max, criterion=criterion, alpha.max.k=alpha.max.k, sigma.sq.max.k=sigma.sq.max.k, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, start.vals=start.vals, mserr=mserr, ub=ub, lb=lb, opts=opts)
    
    cat("Finished. Summarizing", "\n")
    
    ###### NEED TO ADD SUMMARIZING STEP ######
    # Output will have two objects -- the shift locations and the "best" fit, which includes index matrix and criterion value used.
    ######
    
    obj <- list(...)
    return(obj)
}


GetShiftModel <- function(phy, data, nmax, criterion=c("AIC", "AICc", "BIC", "mBIC"), alpha.max.k, sigma.sq.max.k, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, shift.point=0.5, start.vals, mserr="none", ub=20, lb=-21, opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5)){
    
    phy <- reorder(phy, "pruningwise")
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- n + 1L
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]

    curmodel = list()
    flag = 0

    promodel = curmodel
    pro <- NULL
    pro$fit.object <- NULL
    pro$criterion <- Inf

    while (flag==0) {
        for (i in 1:N) {
            if (sum(which(curmodel==i))>0) {
                pos = which(curmodel==i)
                tempmodel = curmodel[-pos]

                ###### Make sure the right stuff is passed here ######
                temp <- OptimizeDredgeLikelihood(tempmodel, phy=phy, data=data, criterion=criterion, alpha.max.k=alpha.max.k, sigma.sq.max.k=sigma.sq.max.k, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, start.vals=start.vals, mserr=mserr, ub=ub, lb=lb, opts=opts)
                #####################################################

                if (temp$criterion < pro$criterion) {
                    promodel = tempmodel
                    pro = temp
                    flag = flag + 1
                }
            }
            
            if (sum(which(curmodel==i))==0) {
                if (length(curmodel) < nmax) {
                    tempmodel = curmodel
                    tempmodel[[length(tempmodel)+1]] = i

                    ###### Make sure the right stuff is passed here ######
                    temp <- OptimizeDredgeLikelihood(tempmodel, phy=phy, data=data, criterion=criterion, alpha.max.k=alpha.max.k, sigma.sq.max.k=sigma.sq.max.k, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, start.vals=start.vals, mserr=mserr, ub=ub, lb=lb, opts=opts)
                    #####################################################
                    print(temp)
                    if (temp$criterion < pro$criterion) {
                        promodel = tempmodel
                        pro = temp
                        flag = flag + 1
                    }
                }
                if (length(curmodel)>=1)
                for (j in 1:length(curmodel))
                if ((anc[i]==des[curmodel[[j]]])||(des[i]==anc[curmodel[[j]]])||(anc[i]==anc[curmodel[[j]]])) {
                    tempmodel = curmodel
                    tempmodel[j] = i

                    ###### Make sure the right stuff is passed here ######
                    temp <- OptimizeDredgeLikelihood(tempmodel, phy=phy, data=data, criterion=criterion, alpha.max.k=alpha.max.k, sigma.sq.max.k=sigma.sq.max.k, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, start.vals=start.vals, mserr=mserr, ub=ub, lb=lb, opts=opts)
                    #####################################################
                    
                    if (temp$criterion < pro$criterion) {
                        promodel = tempmodel
                        pro = temp
                        flag = flag + 1
                    }
                }
            }
        }
        if (flag > 0) flag = 0 else flag = 1
        curmodel <- promodel
        current <- pro
    }
    maps.and.pars <- list(shiftmodel=curmodel, model.fit=current)
    return(maps.and.pars)
}

## TEST ##
#dat <- data.frame(taxon = flowerTree$tip.label, trait=flowerSize$log_transformed_size)
#phy <- flowerTree
#GetShiftModel(phy, dat, nmax=1, criterion="mBIC", alpha.max.k=1, sigma.sq.max.k=1, start.vals=c(0.005824766, 0.0119683))


GetShiftMap <- function(curmodel, phy, data){
    phy = reorder(phy, "pruningwise")
    regimes <- length(curmodel)
    data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]
    new.painting <- rep(1, Nnode(phy) + Ntip(phy))
    for(index in 1:regimes){
        nodes.to.paint <- c( phy$edge[curmodel[[index]],2], getDescendants(phy, phy$edge[curmodel[[index]],2], curr=NULL))
        new.painting[nodes.to.paint] <- index + 1
    }
    TIPS <- 1:Ntip(phy)
    new.data <- data.frame(taxon=phy$tip.label, regime=new.painting[TIPS], trait=data.new[,2])
    phy$node.label <- new.painting[-TIPS]
    new.map <- list(phy=phy, data=new.data)
    return(new.map)
}

### TEST ###
#dat <- cbind(flowerTree$tip.label, flowerSize$log_transformed_size)
#curmodel <- list()
#curmodel[[1]] <- 45
#ll <- GetShiftMap(curmodel, flowerTree, dat)
#plot(ll$phy)
#nodelabels(ll$phy$node.label)
#The model number is the ROW of the pruningwise edge matrix[,2].


DredgeCombinations <- function(shifts, alpha.max.k, sigma.sq.max.k, start.vals){
    
    if(shifts == 1){
        if(alpha.max.k > 1){
            alpha <- expand.grid(1,1:2)
        }else{
            alpha <- expand.grid(1, 1)
        }
        if(sigma.sq.max.k > 1){
            sigma.sq <- expand.grid(1, 1:2)
        }else{
            sigma.sq <- expand.grid(1, 1)
        }
    }
    
    if(shifts == 2){
        if(alpha.max.k > 1){
            alpha <- expand.grid(1, 1:2, 1:3)[-5,]
        }else{
            alpha <- expand.grid(1, 1)
        }
        if(sigma.sq.max.k > 1){
            sigma.sq <- expand.grid(1, 1:2, 1:3)[-5,]
        }else{
            sigma.sq <- expand.grid(1, 1)
        }
    }

    if(shifts == 3){
        if(alpha.max.k > 1){
            alpha <- expand.grid(1, 1:2, 1:3, 1:4)[-c(5,13,17,19:23),]
        }else{
            alpha <- expand.grid(1, 1)
        }
        if(sigma.sq.max.k > 1){
            sigma.sq <- expand.grid(1, 1:2, 1:3, 1:4)[-c(5,13,17,19:23),]
        }else{
            sigma.sq <- expand.grid(1, 1)
        }
    }
    combos <- list(alpha.par.map <- alpha, sigma.sq.par.map <- sigma.sq)
    return(combos)
}

### TEST ###
#DredgeCombinations(shifts=3, alpha.max.k=3, sigma.sq.max.k=3)
#The model number is the ROW of the pruningwise edge matrix[,2].


GetLikelihood <- function(p, phy, data, simmap.tree, root.age, scaleHeight, root.station, shift.point, index.mat=index.mat, mserr=mserr){
    p.new <- exp(p)
    Rate.mat <- index.mat
    Rate.mat[] <- c(p.new, 1e-10)[index.mat]
    tmp <- NA
    try(tmp <- OUwie::OUwie.fixed(phy=phy, data=data, model=c("OUMVA"), simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, alpha=Rate.mat[1,], sigma.sq=Rate.mat[2,], theta=NULL, mserr=mserr, check.identify=FALSE, quiet=TRUE)$loglik, silent=TRUE)
    if(!is.finite(tmp)){
        return(10000000)
    }
    if(is.na(tmp)){
        return(10000000)
    }else{
        loglik <- -1*tmp
    }
    return(loglik)
}


OptimizeDredgeLikelihood <- function(curmodel, phy, data, criterion=c("AIC", "AICc", "BIC", "mBIC"), alpha.max.k, sigma.sq.max.k, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, shift.point=0.5, start.vals, mserr="none", ub=20, lb=-21, opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5)){
    
    ntips <- Ntip(phy)
    if(length(curmodel)==0){
        best.fit <- NULL
        best.fit$fit.object <- start.vals
        if(criterion == "AIC"){
            best.fit$criterion <- start.vals$AIC
        }
        if(criterion == "AICc"){
            best.fit$criterion <- start.vals$AICc
        }
        if(criterion == "BIC"){
            best.fit$criterion <- start.vals$BIC
        }
        if(criterion == "mBIC"){
            best.fit$criterion <- start.vals$logl + log(ntips)
        }

        return(best.fit)
    }
    
    ####GENERATE TREE PAINTING AND DATA SET####
    mapping.tree.data <- GetShiftMap(curmodel, phy, data)
    shifts <- k <- length(curmodel)
    
    ####GET MODEL COMBINATIONS###
    dredge.combos <- DredgeCombinations(shifts=shifts, alpha.max.k=alpha.max.k, sigma.sq.max.k=sigma.sq.max.k)
    ###########################################

    check.identify <- OUwie::check.identify(phy=mapping.tree.data$phy, data=mapping.tree.data$data, simmap.tree=FALSE, get.penalty=TRUE, quiet=TRUE)
    if(check.identify[1] == 0){
        best.fit <- NULL
        best.fit$fit.object <- NULL
        best.fit$criterion <- Inf
        return(best.fit)
    }
    
    index.mat <- def.set.pars <- matrix(0, 2, k+1)
    def.set.pars <- c(rep(start.vals$solution[1,1], k+1), rep(start.vals$solution[2,1], k+1))
    
    current.best.score <- Inf
    current.fit <- NULL
    
    for(alpha.index in 1:dim(dredge.combos$alpha.par.map)[1]){
        for(sigma.index in 1:dim(dredge.combos$sigma.sq.par.map)[1]){
            index.mat[1,] <- as.numeric(dredge.combos$alpha.par.map[alpha.index,])
            max.par <- max(index.mat[1,])
            index.mat[2,] <- as.numeric(dredge.combos$sigma.sq.par.map[sigma.index,]) + max.par
            np <- max(index.mat)
            pars <- c(index.mat[1,], index.mat[2,])
            param.count <- np + (k+1)
            if(root.station == FALSE){
                param.count <- param.count + 1
            }
            
            np.sequence <- 1:np
            ip <- numeric(np)
            upper <- numeric(np)
            for(i in np.sequence){
                ip[i] <- def.set.pars[which(pars == np.sequence[i])[1]]
            }
            
            lower <- rep(lb, length(ip))
            upper <- rep(ub, length(ip))
            
            if(mserr=="est"){
                ip <- c(ip,start.vals$mserr)
                lower <- c(lower,0)
                upper <- c(upper,ub)
            }

            #optimize likelihood
            out <- nloptr(x0=log(ip), eval_f=GetLikelihood, phy=mapping.tree.data$phy, data=mapping.tree.data$data, simmap.tree=FALSE, root.age=NULL, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, index.mat=index.mat, mserr=mserr, lb=lower, ub=upper, opts=opts)
            loglik <- out$objective
            
            if(criterion == "AIC"){
                score <- (2*loglik) + (2*param.count)
            }
            if(criterion == "AICc"){
                score <- (2*loglik) + (2*param.count*(ntips/(ntips-param.count-1)))
            }
            if(criterion == "BIC"){
                score <- (2*loglik) + (log(ntips) * param.count)
            }
            if(criterion == "mBIC"){
                score <- (2*loglik) + check.identify[2]
            }
            if(score < current.best.score){
                current.best.score <- score
                current.fit <- out
            }
        }
    }
    best.fit <- list(fit.object=current.fit, criterion=current.best.score, index.mat=index.mat)
    return(best.fit)
}

### TEST ###
#dat <- data.frame(taxon = flowerTree$tip.label, trait=flowerSize$log_transformed_size)
#curmodel <- list()
#curmodel[[1]] <- 45
#OptimizeDredgeLikelihood(curmodel, phy=flowerTree, data=dat, alpha.max.k=3, sigma.sq.max.k=3, root.station=FALSE, start.vals=c(0.005824766, 0.0119683), mserr="none", opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5))































#this creates a vector. The third element in this vector corresponds to the number appearing in the
#   second column (descendant) of the edge matrix in the phy object for that edge.
#   So if you want to know the ancestor of the i-th element in the rgenoud.individual,
#   the focal node number is phy$edge[(get.mapping(phy,rgenoud.individual))[i],2]
#   and the parent node number is phy$edge[(get.mapping(phy,rgenoud.individual))[i],1]
#   the state of the parent node in the rgenoud.individual is then
#   rgenoud.individual[ which(get.mapping(phy,rgenoud.individual) == phy$edge[(get.mapping(phy,rgenoud.individual))[i],1]) ]
get.mapping <- function(rgenoud.individual,phy) {
    mapping <- match(phy$edge[,2],1:length(rgenoud.individual))
    mapping <- append(mapping,which(!(sequence(max(mapping)) %in% mapping))) #add on the root taxon
    return(mapping)
}


get.final.label<-function(i,rgenoud.individual,phy) {
    if (rgenoud.individual[i]!=0) {
        return(rgenoud.individual[i])
    }
    indexOffset<-(length(rgenoud.individual)/3) * floor(i/(length(rgenoud.individual)/3))
    nodeCount <- Nnode(phy,internal.only=FALSE)
    mapping <- get.mapping(rgenoud.individual,phy)
    while (rgenoud.individual[i] == 0 ) {
        indexMapping <- i %% nodeCount
        current.node <- mapping[ indexMapping ]
        parent.node <- phy$edge[ which(phy$edge[,2] == current.node), 1]
        i <- which( mapping == parent.node ) + indexOffset
        if (length(i) ==0 ) {
            return ( 1) #cannot find the parent node because we ARE at the root
        }
    }
    return( rgenoud.individual[i] )
}



as.full.regime<-function(rgenoud.individual,phy,alphaZero=FALSE) {
    rgenoud.individual <- sapply(X=sequence(length(rgenoud.individual)),FUN=get.final.label,rgenoud.individual=rgenoud.individual,phy=phy)
    if(alphaZero) {
        rgenoud.individual[which(rgenoud.individual==(-1))] <- 0
    }
    return(rgenoud.individual)
}


GetTrees <- function(x){
    obj<-NULL
    x$rgenoud.individual[x$rgenoud.individual==(-1)] <- 0
    tot <- length(x$rgenoud.individual)/3
    nb.tip <- Ntip(x$phy)
    nb.node <- Nnode(x$phy)
    
    ##Gets node labels for theta
    theta.regimes<-x$rgenoud.individual[1:tot]
    n.labels <- numeric(Nnode(x$phy))
    n.labels[1]<-theta.regimes[tot]
    n.labels[2:Nnode(x$phy)] <- theta.regimes[match((2+Ntip(x$phy)):(Nedge(x$phy)+1),x$phy$edge[,2])]
    x$phy$node.label <- n.labels
    obj$theta <- x$phy
    
    ##Gets node labels for sigma
    sigma.regimes <- x$rgenoud.individual[(tot+1):(2*tot)]
    n.labels <- numeric(Nnode(x$phy))
    n.labels[1] <- sigma.regimes[tot]
    n.labels[2:Nnode(x$phy)] <- sigma.regimes[match((2+Ntip(x$phy)):(Nedge(x$phy)+1),x$phy$edge[,2])]
    x$phy$node.label <- n.labels
    obj$sigma <- x$phy
    
    ##Gets node labels for alpha
    alpha.regimes <- x$rgenoud.individual[(2*tot+1):length(x$rgenoud.individual)]
    n.labels <- numeric(Nnode(x$phy))
    n.labels[1] <- alpha.regimes[tot]
    n.labels[2:Nnode(x$phy)] <- alpha.regimes[match((2+Ntip(x$phy)):(Nedge(x$phy)+1),x$phy$edge[,2])]
    x$phy$node.label <- n.labels
    obj$alpha <- x$phy
    
    obj
}


print.ouwie.dredge.result <- function(x, ...) {
    K<-sum(unlist(param.count(x$rgenoud.individual,x$phy)))
    n<-Ntip(x$phy)
    output<-c(x$loglik,x$AIC,x$AICc,n,dim(x$regime.mat)[1],K,c(unlist(param.count(x$rgenoud.individual,x$phy))))
    names(output) <- c("negLnL","AIC","AICc", "mBIC", "ntax", "n_regimes", "K_all", "K_theta", "K_sigma", "K_alpha")
    print(output)
    
    cat("\nRegime matrix\n")
    colnames(x$regime.mat)<-c("theta","sigma","alpha")
    rownames(x$regime.mat)<-c(paste("regime_",sequence(dim(x$regime.mat)[1]),sep=""))
    print(x$regime.mat)
    
    regime.mat.params<-x$regime.mat*0
    p.index<-1
    for(k in sequence(param.count(x$rgenoud.individual,x$phy)[[1]])) {
        regime.mat.params[which(x$regime.mat[,1]==k),1]<-x$solution[p.index]
        p.index<-p.index+1
    }
    for(k in sequence(param.count(x$rgenoud.individual,x$phy)[[2]])) {
        regime.mat.params[which(x$regime.mat[,2]==k),2]<-x$solution[p.index]
        p.index<-p.index+1
    }
    for(k in sequence(param.count(x$rgenoud.individual,x$phy)[[3]])) {
        regime.mat.params[which(x$regime.mat[,3]==k),3]<-x$solution[p.index]
        p.index<-p.index+1
    }
    cat("\nRate matrix\n")
    print(regime.mat.params)
    cat("\n")
    
}

##   Regime vs. Rate        ##
##   gradient of color      ##

plot.ouwie.dredge.result <- function(x, type="regime", col.pal=c("Set1"), ...) {
    #    x$rgenoud.individual<-as.full.regime(x$rgenoud.individual,x$phy)
    par(mfcol=c(1,3))
    tot<-length(x$rgenoud.individual)/3
    
    if(type=="regime"){
        
        regimes<-x$rgenoud.individual
        regime.lvls<-length(levels(as.factor(regimes)))
        if(regime.lvls<3){
            regime.lvls = 3
        }
        #
        if(length(col.pal>1)){
            co<-brewer.pal(regime.lvls, col.pal)
        }
        #If not the
        else{
            co<-col.pal
        }
        ##Plot thetas
        nb.tip <- Ntip(x$phy)
        nb.node <- Nnode(x$phy)
        comp <- numeric(Nedge(x$phy))
        comp[match(1:(Ntip(x$phy)), x$phy$edge[,2])] <- as.factor(regimes[1:nb.tip])
        comp[match((2+Ntip(x$phy)):(Nedge(x$phy)+1), x$phy$edge[,2])] <- as.factor(regimes[(nb.tip+1):(nb.tip+nb.node)][-1])
        plot.phylo(x$phy,edge.color=co[comp], ...)
        title(main="theta")
        
        ##Plot sigmas
        nb.tip <- Ntip(x$phy)
        nb.node <- Nnode(x$phy)
        comp <- numeric(Nedge(x$phy))
        comp[match(1:(Ntip(x$phy)), x$phy$edge[,2])] <- as.factor(regimes[(tot+1):(tot+nb.tip)])
        comp[match((2+Ntip(x$phy)):(Nedge(x$phy)+1), x$phy$edge[,2])] <- as.factor(regimes[(2*tot+1):(2*(tot))][-1])
        plot.phylo(x$phy,edge.color=co[comp], ...)
        title(main="sigma")
        
        ##Plot alphas
        nb.tip <- Ntip(x$phy)
        nb.node <- Nnode(x$phy)
        comp <- numeric(Nedge(x$phy))
        comp[match(1:(Ntip(x$phy)), x$phy$edge[,2])] <- as.factor(regimes[(2*tot+1):(2*tot+nb.tip)])
        comp[match((2+Ntip(x$phy)):(Nedge(x$phy)+1), x$phy$edge[,2])] <- as.factor(regimes[(3*tot+1):length(x$rgenoud.individual)][-1])
        plot.phylo(x$phy,edge.color=co[comp], ...)
        title(main="alpha")
    }
}





