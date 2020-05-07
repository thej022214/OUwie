#library(phytools)
#library(OUwie)
#need add getDescendants to the NAMESPACE
#library(RColorBrewer)


OUwie.dredge <- function(phy, data, criterion=c("AIC", "AICc", "BIC", "mBIC"), shift.max=3, sigma.sq.max.k=3, alpha.max.k=3, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, shift.point=0.5, mserr="none", opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5)) {
    
    ### ADD WARNINGS ###
    ## Number of alpha or sigma pars cannot exceed 1+max shifts
    ##
    
    #Coerce the data so that it will run in OUwie -- using values of OU1 as the starting points:
    cat("Initializing...","\n")
    
    data2 <- data.frame(taxon=as.character(data[,1]), regime=rep(1, length(data[,1])), trait=data[,2], stringsAsFactors=FALSE)
    #phy$node.label <- sample(c(1:2),phy$Nnode, replace=TRUE)
    start.vals <- OUwie(phy, data2, model=c("OU1"), quiet=TRUE, root.station=TRUE, scaleHeight=scaleHeight, mserr=mserr, check.identify=FALSE)
    cat("Begin optimization routine -- Starting values:", c(start.vals$solution[1,1], start.vals$solution[2,1]), "\n")
    phy$node.label <- NULL
    find.shifts <- GetShiftModel(phy=phy, data=data, nmax=shift.max, criterion=criterion, alpha.max.k=alpha.max.k, sigma.sq.max.k=sigma.sq.max.k, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, start.vals=start.vals, mserr=mserr, opts=opts)

    cat("Finished. Summarizing", "\n")

    #Step 1: Paint the tree up and make the proper data.
    mapped.phy <- find.shifts$model.fit$model.phy
    if(dim(find.shifts$model.fit$model.data)[2] == 2){
        mapped.data <- cbind(phy$tip.label, find.shifts$model.fit$model.data)
    }else{
        mapped.data <- find.shifts$model.fit$model.data
    }
    
    #Step 2: Get thetas:
    mapped.thetas <- GetThetas(p=log(find.shifts$model.fit$fit.object$solution), phy=mapped.phy, data=mapped.data, simmap.tree=FALSE, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, index.mat=find.shifts$model.fit$index.mat, mserr=mserr)
    regime.weights <- mapped.thetas$regime.weights

    #Step 3: Now summarize everything:
    solution <- mapped.thetas$solution
    rownames(solution) <- rownames(find.shifts$model.fit$index.mat) <- c("alpha", "sigma.sq")
    tot.states <- factor(c(mapped.phy$node.label, as.character(mapped.data[,2])))
    colnames(solution) <- levels(tot.states)

    if(mserr=="est"){
        mserr.est <- find.shifts$solution[length(find.shifts$solution)]
        param.count <- find.shifts$param.count
    }
    else{
        mserr.est <- NULL
    }
    
    obj <- list(loglik = mapped.thetas$loglik, criterion=criterion, criterion.score=find.shifts$model.fit$criterion, shift.model=find.shifts$shiftmodel, solution=solution, mserr.est=mserr.est, theta=mapped.thetas$theta, tot.states=tot.states, index.mat=find.shifts$model.fit$fit.object$index.mat, simmap.tree=FALSE, root.age=root.age, scaleHeight=scaleHeight, shift.point=shift.point, opts=opts, data=mapped.data, phy=mapped.phy, root.station=root.station, starting.vals=start.vals$solution, regime.weights=regime.weights)
    class(obj) <- "OUwie.dredge"
    return(obj)
}


print.OUwie.dredge <- function(x, ...) {
    
    n <- Ntip(x$phy)

    output <- data.frame(x$loglik, x$criterion, x$criterion.score, n, length(x$shift.point), length(unique(x$solution[1,])), length(unique(x$solution[2,])), row.names="")
    names(output) <- c("negLnL", "criterion", "score", "ntax", "n_regimes", "k_alpha", "k_sigma")
    
    cat("\nFit\n")
    print(output)
    cat("\n")

    cat("\nRegime matrix\n")
    last.column.to.remove <- dim(x$regime.weights)[2]
    print(x$regime.weights[1:3, -last.column.to.remove])
    
    cat("\n")
    
}


GetShiftModel <- function(phy, data, nmax, criterion=c("AIC", "AICc", "BIC", "mBIC"), alpha.max.k, sigma.sq.max.k, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, shift.point=0.5, start.vals, mserr="none", opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5)){
    
    phy <- reorder(phy, "pruningwise")
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
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
                temp <- OptimizeDredgeLikelihood(curmodel=tempmodel, phy=phy, data=data, criterion=criterion, alpha.max.k=alpha.max.k, sigma.sq.max.k=sigma.sq.max.k, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, start.vals=start.vals, mserr=mserr, opts=opts)
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
                    temp <- OptimizeDredgeLikelihood(curmodel=tempmodel, phy=phy, data=data, criterion=criterion, alpha.max.k=alpha.max.k, sigma.sq.max.k=sigma.sq.max.k, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, start.vals=start.vals, mserr=mserr, opts=opts)
                    #####################################################

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
                    temp <- OptimizeDredgeLikelihood(curmodel=tempmodel, phy=phy, data=data, criterion=criterion, alpha.max.k=alpha.max.k, sigma.sq.max.k=sigma.sq.max.k, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, start.vals=start.vals, mserr=mserr, opts=opts)
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
#library(phylolm)
#data(flowerSize)
#data(flowerTree)
#dat <- data.frame(taxon = flowerTree$tip.label, trait=flowerSize$log_transformed_size)
#phy <- flowerTree
#start.values <- start.vals
#GetShiftModel(phy, dat, nmax=1, criterion="mBIC", alpha.max.k=1, sigma.sq.max.k=1, start.vals=start.values, root.station=TRUE)


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
            alpha <- expand.grid(1, 1, 1)
        }
        if(sigma.sq.max.k > 1){
            sigma.sq <- expand.grid(1, 1:2, 1:3)[-5,]
        }else{
            sigma.sq <- expand.grid(1, 1, 1)
        }
    }

    if(shifts == 3){
        if(alpha.max.k > 1){
            alpha <- expand.grid(1, 1:2, 1:3, 1:4)[-c(5,13,17,19:23),]
        }else{
            alpha <- expand.grid(1, 1, 1, 1)
        }
        if(sigma.sq.max.k > 1){
            sigma.sq <- expand.grid(1, 1:2, 1:3, 1:4)[-c(5,13,17,19:23),]
        }else{
            sigma.sq <- expand.grid(1, 1, 1, 1)
        }
    }
    combos <- list(alpha.par.map=alpha, sigma.sq.par.map=sigma.sq)
    return(combos)
}

### TEST ###
#DredgeCombinations(shifts=3, alpha.max.k=3, sigma.sq.max.k=3)
#The model number is the ROW of the pruningwise edge matrix[,2].
GetThetas <- function(p, phy, data, simmap.tree, root.age, scaleHeight, root.station, shift.point, index.mat=index.mat, mserr=mserr){
    p.new <- exp(p)
    Rate.mat <- index.mat
    Rate.mat[] <- c(p.new, 1e-10)[index.mat]
    tmp <- NA
    try(tmp <- OUwie.fixed(phy=phy, data=data, model=c("OUM"), simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, alpha=Rate.mat[1,], sigma.sq=Rate.mat[2,], theta=NULL, mserr=mserr, check.identify=FALSE, quiet=TRUE), silent=TRUE)
    if(!is.finite(tmp[[1]])){
        return(NULL)
    }
    if(is.na(tmp[[1]])){
        return(NULL)
    }else{
        return(tmp)
    }
}


GetLikelihood <- function(p, phy, data, simmap.tree, root.age, scaleHeight, root.station, shift.point, index.mat=index.mat, mserr=mserr){
    p.new <- exp(p)
    Rate.mat <- index.mat
    Rate.mat[] <- c(p.new, 1e-10)[index.mat]
    tmp <- NA
    try(tmp <- OUwie.fixed(phy=phy, data=data, model=c("OUM"), simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, alpha=Rate.mat[1,], sigma.sq=Rate.mat[2,], theta=NULL, mserr=mserr, check.identify=FALSE, quiet=TRUE)$loglik, silent=TRUE)
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


OptimizeDredgeLikelihood <- function(curmodel, phy, data, criterion=c("AIC", "AICc", "BIC", "mBIC"), alpha.max.k, sigma.sq.max.k, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, shift.point=0.5, start.vals, mserr="none", opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5)){
   
    ub=20
    lb=-21
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
            current.best.score <- (-2 * start.vals$loglik) + log(ntips)
            best.fit <- list(fit.object=start.vals, criterion=current.best.score, index.mat=start.vals$index.mat, param.count=start.vals$param.count, model.phy=start.vals$phy, model.data=start.vals$data)
        }
        return(best.fit)
    }
    
    ####GENERATE TREE PAINTING AND DATA SET####
    mapping.tree.data <- GetShiftMap(curmodel, phy, data)
    shifts <- k <- length(curmodel)
    ####GET MODEL COMBINATIONS###
    dredge.combos <- DredgeCombinations(shifts=shifts, alpha.max.k=alpha.max.k, sigma.sq.max.k=sigma.sq.max.k)
    ###########################################

    check.identify <- check.identify.dredge(phy=mapping.tree.data$phy, data=mapping.tree.data$data, simmap.tree=FALSE, get.penalty=TRUE, quiet=TRUE)
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
            param.count <- np
            if(root.station == FALSE){
                param.count <- param.count
            }
            
            np.sequence <- 1:np
            ip <- numeric(np)
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
            out <- nloptr(x0=log(ip), eval_f=GetLikelihood, phy=mapping.tree.data$phy, data=mapping.tree.data$data, simmap.tree=FALSE, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, shift.point=shift.point, index.mat=index.mat, mserr=mserr, lb=lower, ub=upper, opts=opts)
            out$solution <- exp(out$solution)
            loglik <- out$objective
            if(criterion == "AIC"){
                score <- (2*loglik) + (2*(shifts+param.count))
            }
            if(criterion == "sAICc"){
                #What surface does
                score <- (2*loglik) + (2*(2 * shifts + param.count))
            }
            if(criterion == "BIC"){
                score <- (2*loglik) + (log(ntips)*(shifts + param.count))
            }
            if(criterion == "mBIC"){
                #need to finish penalty calculation --
                penalty <- (param.count * log(ntips)) + check.identify[2]
                score <- (2*loglik) + penalty
            }
            if(score < current.best.score){
                current.best.score <- score
                current.fit <- out
                current.index.mat <- index.mat
                current.param.count <- param.count
                current.data <- mapping.tree.data$data
                current.phy <- mapping.tree.data$phy
            }
        }
    }
    best.fit <- list(fit.object=current.fit, criterion=current.best.score, index.mat=current.index.mat, param.count=current.param.count, model.phy=current.phy, model.data=current.data)
    return(best.fit)
}

### TEST ###
#library(phylolm)
#dat <- data.frame(taxon = flowerTree$tip.label, trait=flowerSize$log_transformed_size)
#curmodel <- list()
#curmodel[[1]] <- 45
#OptimizeDredgeLikelihood(curmodel, phy=flowerTree, data=dat, alpha.max.k=3, sigma.sq.max.k=3, root.station=FALSE, start.vals=c(0.005824766, 0.0119683), mserr="none", opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5))

GetParameterPainting <- function(phy, data, rates, k.pars){
    full.regime.map <- c(data[,2], phy$node.label)
    regimes <- length(unique(c(data[,2], phy$node.label)))
    painting <- rep(1, Nnode(phy) + Ntip(phy))
    if(k.pars == 1){
        painting <- rep(1, Nnode(phy) + Ntip(phy))
    }else{
        for(index in 2:regimes){
            if(rates[index] != rates[c(index-1)]){
                painting[full.regime.map==index] <- index
            }
        }
    }
    return(painting)
}


plot.OUwie.dredge <- function(x, col.pal=c("Set1"), ...) {
    
    truncated.regime.weight.mat <- x$regime.weights[,-c(dim(x$regime.weights)[2])]
    alpha.regimes <- truncated.regime.weight.mat[1,]
    sigma.regimes <- truncated.regime.weight.mat[2,]
    theta.regimes <- truncated.regime.weight.mat[3,]
     
    par(mfcol=c(1,3))
    
    regime.lvls <- dim(truncated.regime.weight.mat)[2]
    if(regime.lvls<3){
        regime.lvls = 3
    }
    if(length(col.pal>1)){
        co <- brewer.pal(regime.lvls, col.pal)
    }else{
        co <- col.pal
    }
    
    ##Plot thetas
    regimes <- GetParameterPainting(phy=x$phy, data=x$data, rates=truncated.regime.weight.mat[3,], k.pars=length(unique(truncated.regime.weight.mat[3,])))
    nb.tip <- Ntip(x$phy)
    nb.node <- Nnode(x$phy)
    comp <- numeric(Nedge(x$phy))
    comp[match(1:(Ntip(x$phy)), x$phy$edge[,2])] <- as.factor(regimes[1:nb.tip])
    comp[match((2+Ntip(x$phy)):(Nedge(x$phy)+1), x$phy$edge[,2])] <- as.factor(regimes[(nb.tip+1):(nb.tip+nb.node)][-1])
    plot.phylo(x$phy, edge.color=co[comp], ...)
    title(main=expression(theta))
    
    ##Plot sigmas
    regimes <- GetParameterPainting(phy=x$phy, data=x$data, rates=truncated.regime.weight.mat[2,], k.pars=length(unique(truncated.regime.weight.mat[2,])))
    nb.tip <- Ntip(x$phy)
    nb.node <- Nnode(x$phy)
    comp <- numeric(Nedge(x$phy))
    comp[match(1:(Ntip(x$phy)), x$phy$edge[,2])] <- as.factor(regimes[1:nb.tip])
    comp[match((2+Ntip(x$phy)):(Nedge(x$phy)+1), x$phy$edge[,2])] <- as.factor(regimes[(nb.tip+1):(nb.tip+nb.node)][-1])
    plot.phylo(x$phy, edge.color=co[comp], ...)
    title(main=expression(sigma^2))
    
    ##Plot alphas
    regimes <- GetParameterPainting(phy=x$phy, data=x$data, rates=truncated.regime.weight.mat[1,], k.pars=length(unique(truncated.regime.weight.mat[1,])))
    nb.tip <- Ntip(x$phy)
    nb.node <- Nnode(x$phy)
    comp <- numeric(Nedge(x$phy))
    comp[match(1:(Ntip(x$phy)), x$phy$edge[,2])] <- as.factor(regimes[1:nb.tip])
    comp[match((2+Ntip(x$phy)):(Nedge(x$phy)+1), x$phy$edge[,2])] <- as.factor(regimes[(nb.tip+1):(nb.tip+nb.node)][-1])
    plot.phylo(x$phy, edge.color=co[comp], ...)
    title(main=expression(alpha))
    
}





