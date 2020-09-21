
######################################################################################################################################
######################################################################################################################################
### Contour plots OUwie analyses
######################################################################################################################################
######################################################################################################################################

#written by Jeremy M. Beaulieu

contourSearchPoints <- function(variables, lower, upper, nreps){
    #Creates a latin square design of the parameter space:
    X <- randomLHS(nreps, variables)
    param.points <- matrix(0, nrow=nreps, ncol=variables)
    param.points[,1] <- qunif(X[,1], lower[1], upper[1])
    param.points[,2] <- qunif(X[,2], lower[2], upper[2])
    return(param.points)
}


OUwie.semifixed <- function(p, phy, data, num.regimes, index.vector, fixed.pars, simmap.tree, scaleHeight, root.station, get.root.theta, shift.point, algorithm){
    
    new.p <- exp(p)
    
    rates <- index.vector
    rates[] <- c(new.p, 42)[index.vector]
    rates[rates==42] <- fixed.pars
    
    alpha <- rates[1:num.regimes]
    rates <- rates[-c(1:num.regimes)]
    sigma.sq <- rates[1:num.regimes]
    rates <- rates[-c(1:num.regimes)]
    theta <- rates
    pp <- NA
    try(pp <- OUwie.fixed(phy=phy, data=data, model="OUMVA", simmap.tree=simmap.tree, scaleHeight=scaleHeight, root.station=root.station, get.root.theta=get.root.theta, shift.point=shift.point, alpha=alpha, sigma.sq=sigma.sq, theta=theta, check.identify=FALSE, algorithm=algorithm, quiet=TRUE), silent=TRUE)
    if(is.na(pp[1])){
        return(1000000)
    }else{
        return(-pp$loglik)
    }
}


contourSearchOUwie <- function(phy, data, num.regimes, param.points, index.vector, par.vector, simmap.tree, scaleHeight, root.station, get.root.theta, shift.point, algorithm, opts, n.cores) {
    
    nreps <- dim(param.points)[1]
    max.pars <- max(index.vector)
    init.vals <- unique(par.vector[index.vector<max.pars])
    lower <- rep(-21, length(init.vals))
    upper <- rep(21, length(init.vals))
    
    if(is.null(n.cores)){
        res <- matrix(,nreps,3)
        for(nrep.index in 1:nreps){
            print(nrep.index)
            fixed.pars <- c(param.points[nrep.index,1], param.points[nrep.index,2])
            opts <- opts
            out <- nloptr(x0=log(init.vals), eval_f=OUwie.semifixed, opts=opts, lb=lower, ub=upper, phy=phy, data=data, num.regimes=num.regimes, index.vector=index.vector, fixed.pars=fixed.pars, simmap.tree=simmap.tree, scaleHeight=scaleHeight, root.station=root.station, get.root.theta=get.root.theta, shift.point=shift.point, algorithm=algorithm)
            res[nrep.index,] <- c(-out$objective, fixed.pars[1], fixed.pars[2])
        }
    }else{
        PointEval <- function(nrep.index){
            fixed.pars <- c(param.points[nrep.index,1], param.points[nrep.index,2])
            opts <- opts
            out <- nloptr(x0=log(init.vals), eval_f=OUwie.semifixed, opts=opts, lb=lower, ub=upper, phy=phy, data=data, num.regimes=num.regimes, index.vector=index.vector, fixed.pars=fixed.pars, simmap.tree=simmap.tree, scaleHeight=scaleHeight, root.station=root.station, get.root.theta=get.root.theta, shift.point=shift.point, algorithm=algorithm)
            return(c(-out$objective, fixed.pars[1], fixed.pars[2]))
        }
        res.list <- mclapply(1:nreps, PointEval, mc.cores=n.cores)
        res <- matrix(unlist(res.list), ncol = 3, byrow = TRUE)
    }
    return(res)
}



OUwie.contour <- function(OUwie.obj, focal.params=c("alpha_1", "sigma.sq_1"), focal.params.lower=c(0,0), focal.params.upper=c(5,5), nreps=1000, n.cores=NULL){
    
    new.data <- data.frame(taxon=rownames(OUwie.obj$data), regime=OUwie.obj$data[,1], trait=OUwie.obj$data[,2])
    
    #Step 1: Need to generate a vector of what gets fixed and what gets estimated.
    full.par <- OUwie.obj$regime.weights[1:3,]
    full.par <- full.par[,-dim(full.par)[2]]

    #Always two:
    focal.param.split1 <- strsplit(focal.params[1], "_")[[1]]
    focal.param.split2 <- strsplit(focal.params[2], "_")[[1]]
    
    mle.dat <- c(OUwie.obj$loglik, full.par[focal.param.split1[1], focal.param.split1[2]], full.par[focal.param.split2[1], focal.param.split2[2]])
    full.par[focal.param.split1[1], focal.param.split1[2]] <- 0
    full.par[focal.param.split2[1], focal.param.split2[2]] <- 0
    
    if(OUwie.obj$get.root.theta == TRUE){
        num.regimes <- dim(full.par)[2] - 1
        par.vector <- c(full.par[1,2:dim(full.par)[2]],full.par[2, 2:dim(full.par)[2]], full.par[3, 1:dim(full.par)[2]])
        index.vector <- par.vector
        unique.pars <- unique(index.vector[index.vector>0])
        par.count <- 1
        for(par.index in 1:length(unique.pars)){
            index.vector[which(index.vector %in% unique.pars[par.index])] <- par.count
            par.count <- par.count + 1
        }
        index.vector[index.vector == 0] <- par.count
        names(index.vector) <- NULL
    }else{
        num.regimes <- dim(full.par)[2]
        par.vector <- c(full.par[1,1:dim(full.par)[2]],full.par[2, 1:dim(full.par)[2]], full.par[3, 1:dim(full.par)[2]])
        index.vector <- par.vector
        unique.pars <- unique(index.vector[index.vector>0])
        par.count <- 1
        for(par.index in 1:length(unique.pars)){
            index.vector[which(index.vector %in% unique.pars[par.index])] <- par.count
            par.count <- par.count + 1
        }
        index.vector[index.vector == 0] <- par.count
        names(index.vector) <- NULL
    }
    
    #Step 2: Get parameter points to fix using a latin square design:
    param.points <- contourSearchPoints(variables=2, lower=focal.params.lower, upper=focal.params.upper, nreps=nreps)
    
    #Step 3: Call contourSearchOUwie -- ISSUE. What if theta is less than 0? Will deal with later...
    surface.data <- contourSearchOUwie(phy=OUwie.obj$phy, data=new.data, num.regimes=num.regimes, param.points=param.points, index.vector=index.vector, par.vector=par.vector, simmap.tree=OUwie.obj$simmap.tree, scaleHeight=OUwie.obj$scaleHeight, root.station=OUwie.obj$root.station, get.root.theta=OUwie.obj$get.root.theta, shift.point=OUwie.obj$shift.point, algorithm=OUwie.obj$algorithm, opts=OUwie.obj$opts, n.cores=n.cores)
    
    surface.data <- rbind(mle.dat, surface.data, deparse.level=0)
    #Step 4: Make plot with results
    
    obj <- list(surface.data=surface.data, focal.params=focal.params, focal.params.lower=focal.params.lower, focal.params.upper=focal.params.upper)
    class(obj) <- "OUwie.contour"
    return(obj)
}




plot.OUwie.contour <- function(x, mle.point=NULL, levels=c(0:20*0.1), xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, col=grey.colors(21, start=0, end=1), ...){
    
    if(is.null(xlab)){
        xlab = x$focal.params[1]
    }
    if(is.null(ylab)){
        ylab = x$focal.params[2]
    }
    if(is.null(xlim)){
        xlim = c(x$focal.params.lower[1], x$focal.params.upper[1], 1)
    }
    if(is.null(ylim)){
        ylim = c(x$focal.params.lower[2], x$focal.params.upper[2], 1)
    }

    mydata <- data.frame(x=matrix(x$surface.data[,2],ncol=1),y=matrix(x$surface.data[,3],ncol=1),z=matrix(x$surface.data[,1],ncol=1))
    mydata$z <- (-1)*mydata$z
    mydata$z <- mydata$z-min(mydata$z)
    
    interp.res <- interp(x=mydata$x, y=mydata$y, z=mydata$z, xo=seq(min(mydata$x), max(mydata$x),length = 400), yo=seq(min(mydata$y), max(mydata$y),length = 400), duplicate=FALSE)
    plot(NA,xlab="", ylab="", frame=FALSE, axes=FALSE, xaxs="i", yaxs="i", ylim=ylim[1:2], xlim=xlim[1:2], ...)
    .filled.contour(interp.res$x, interp.res$y, interp.res$z, levels=levels, col=col)
    par(tck=.01)
    axis(2, at = seq(ylim[1], ylim[2], by = ylim[3]), las=1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
    axis(1, at = seq(xlim[1], xlim[2], by = xlim[3]), las=1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
    if(!is.null(mle.point)){
        points(x=x$surface.data[1,2], y=x$surface.data[1,3], pch=19, col=mle.point)
    }
    title(xlab=xlab, line=2.5)
    title(ylab=ylab, line=2)
}




#load("simsOUidentify_74")
#oum$scaleHeight <- TRUE
#surf.dat <- OUwie.liksurface(OUwie.obj=oum, simmap.tree=FALSE, focal.param=c("theta_Root", "theta_1"), focal.param.lower=c(0, 0), focal.param.upper=c(10,10), nreps=50)
#levels <- c(0:20*0.1)


