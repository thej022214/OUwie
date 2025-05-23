##OUwie Simulator##

#written by Jeremy M. Beaulieu

#Simulates the Hansen model of continuous characters evolving under discrete selective
#regimes. The input is a tree of class "phylo" that has the regimes as internal node labels
#and a data file the contains the regime states for each species. The trait file must be in
#the following order: Species names then Regime. The user must specify the parameters values
#for each simulation (i.e. alpha, sigma.sq, theta0, theta).

##The following examples assume 2 selective regimes and different models can be specified:
#single rate Brownian motion BM1: alpha=c(0,0); sigma.sq=c(0.9); theta0=0; theta=c(0,0)
#two rate Brownian motion BMS: alpha=c(0,0); sigma.sq=c(0.45,.9); theta0=0; theta=c(0,0)
#global OU (OU1): alpha=c(0.1,0.1); sigma.sq=c(0.9,0.9); theta0=1; theta=c(1,1)
#normal OU (OUSM): alpha=c(0.1,0.1); sigma.sq=c(0.9,0.9); theta0=0; theta=c(1,2)
#multiple sigmas (OUSMV): alpha=c(0.1,0.1); sigma.sq=c(0.45,0.9);
#multiple alphas (OUSMA): alpha=c(0.5,0.1); sigma.sq=c(0.9,0.9); theta0=0; theta=c(1,2)
#multiple alphas and sigmas (OUSMVA): alpha=c(0.5,0.1); sigma.sq=c(0.45,0.9); theta0=0; theta=c(1,2)

OUwie.sim <- function(phy=NULL, data=NULL, simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE, alpha=NULL, sigma.sq=NULL, theta0=NULL, theta=NULL, tip.fog="none", shift.point=0.5, fitted.object=NULL, get.all=FALSE){
	tip.fog_vector <- NA
    if(!is.null(fitted.object)) {
        if(grepl("BM", fitted.object$model) | grepl("OU1", fitted.object$model)) {
            stop(paste("not implemented yet for ", fitted.object$model))
        }
        if(!is.null(alpha) | !is.null(theta0) | !is.null(theta)) {
            stop("You're passing in parameters to simulate from AND a fitted object to simulate under. You can do one or the other")
        }
        phy <- fitted.object$phy
        data <- cbind(phy$tip.label, fitted.object$data)
        alpha <- fitted.object$solution['alpha',]
        alpha[which(is.na(alpha))] <- 0
        sigma.sq <- fitted.object$solution['sigma.sq',]
        
        if(tip.fog != "none"){
            warning("Tip fog is not yet handled for simulations from fitted.object")
        }
        
        if (fitted.object$root.station == TRUE | fitted.object$root.station==FALSE){
            if (fitted.object$model == "OU1"){
                theta <- matrix(t(fitted.object$theta[1,]), 2, length(levels(fitted.object$tot.states)))[1,]
                theta0 <- theta[phy$node.label[1]]
            }
        }
        if (fitted.object$root.station == TRUE | !grepl("OU", fitted.object$model)){ # BM1 or BMS as well
            if (fitted.object$model != "OU1"){
                theta <- matrix(t(fitted.object$theta), 2, length(levels(fitted.object$tot.states)))[1,]
                theta0 <- theta[phy$node.label[1]]
            }
        }
        if (fitted.object$root.station == FALSE & grepl("OU", fitted.object$model)){
            if (fitted.object$model != "OU1"){
                if(fitted.object$get.root.theta == TRUE){
                    theta.all <- matrix(t(fitted.object$theta), 2, 1:length(levels(fitted.object$tot.states))+1)[1,]
                    theta <- theta.all[2:length(theta.all)]
                    theta0 <- theta.all[1]
                }else{
                    int.states <- factor(phy$node.label)
                    phy.tmp <- phy
                    phy.tmp$node.label <- as.numeric(int.states)
                    theta <- matrix(t(fitted.object$theta), 2, length(levels(fitted.object$tot.states)))[1,]
                    theta0 <- theta[phy.tmp$node.label[1]]
                }
            }
        }
    }

    if(is.null(root.age)){
        if(any(branching.times(phy)<0)){
            stop("Looks like your tree is producing negative branching times. Must input known root age of tree.", .call=FALSE)
        }
    }

    #Makes sure the data is in the same order as the tip labels
	if(simmap.tree == FALSE){
		#This is annoying, but the second column has to be in there twice otherwise, error.
        if(tip.fog == "none"){
            data <- data.frame(data[,2], data[,2], row.names=data[,1])
        }
        if(tip.fog == "known"){
            data <- data.frame(data[,2], data[,3], row.names=data[,1])
			tip.fog_vector <- data[,2] #because we've shifted things over
        }
		if(is.numeric(tip.fog)){
			if(length(tip.fog) == length(phy$tip.label)){
				data <- data.frame(data[,2], tip.fog, row.names=data[,1])
				tip.fog_vector <- tip.fog
			}
			if(length(tip.fog)==1){
				data <- data.frame(data[,2], rep(tip.fog, length(phy$tip.label)), row.names=data[,1])
				tip.fog_vector <- rep(tip.fog, length(phy$tip.label))
			}
			tip.fog <- "known"
		}

		data <- data[phy$tip.label,]

		n <- max(phy$edge[,1])
		ntips <- length(phy$tip.label)

		int.states <- factor(phy$node.label)
		phy$node.label <- as.numeric(int.states)
		tip.states <- factor(data[,1])
		data[,1] <- as.numeric(tip.states)
		tot.states <- factor(c(phy$node.label,as.character(data[,1])))
		k <- length(levels(tot.states))

		regime <- matrix(rep(0,(n-1)*k), n-1, k)

		#Obtain root state and internal node labels
		root.state <- phy$node.label[1]
		int.state <- phy$node.label[-1]

		#New tree matrix to be used for subsetting regimes
		edges <- cbind(c(1:(n-1)),phy$edge,MakeAgeTable(phy, root.age=root.age))
		if(scaleHeight == TRUE){
			edges[,4:5] <- edges[,4:5]/max(MakeAgeTable(phy, root.age=root.age))
            root.age <- 1
		}
        
		edges <- edges[sort.list(edges[,3]),]
		mm <- c(data[,1],int.state)

		regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
		#Generates an indicator matrix from the regime vector
		for (i in 1:length(mm)) {
			regime[i,mm[i]] <- 1
		}
		#Finishes the edges matrix
		edges <- cbind(edges,regime)

		#Resort the edge matrix so that it looks like the original matrix order
		edges <- edges[sort.list(edges[,1]),]

		oldregime <- root.state

		alpha <- alpha
		alpha[alpha==0] <- 1e-10
		sigma <- sqrt(sigma.sq)
		theta  <- theta

		x <- matrix(0, n, 1)
		TIPS <- 1:ntips
		ROOT <- ntips + 1L
		x[ROOT,] <- theta0

		for(i in 1:length(edges[,1])){
			anc <- edges[i,2]
			desc <- edges[i,3]
			oldtime <- edges[i,4]
			newtime <- edges[i,5]
			if(anc%in%edges[,3]){
				start <- which(edges[,3]==anc)
				oldregime <- which(edges[start,6:(k+5)]==1)
			}else{
				#For the root:
				oldregime <- root.state
			}
			newregime=which(edges[i,6:(k+5)]==1)

			if(oldregime==newregime){
				x[edges[i,3],] <- (x[edges[i,2],]*exp(-alpha[oldregime]*(newtime-oldtime))) + (theta[oldregime]*(1-exp(-alpha[oldregime]*(newtime-oldtime)))) + (sigma[oldregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[oldregime]*(newtime-oldtime)))/(2*alpha[oldregime])))
			}else{
                shifttime <- newtime-((newtime-oldtime) * shift.point)
				epoch1 <- (x[edges[i,2],]*exp(-alpha[oldregime]*(shifttime-oldtime))) + (theta[oldregime]*(1-exp(-alpha[oldregime]*(shifttime-oldtime)))) + (sigma[oldregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[oldregime]*(shifttime-oldtime)))/(2*alpha[oldregime])))
				oldtime <- shifttime
				newtime <- newtime
				x[edges[i,3],] <- (epoch1*exp(-alpha[newregime]*(newtime-oldtime))) + (theta[newregime]*(1-exp(-alpha[newregime]*(newtime-oldtime)))) + (sigma[newregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[newregime]*(newtime-oldtime)))/(2*alpha[newregime])))
			}
		}

        if(get.all == TRUE){
            sim.dat <- matrix(,length(x),3)
            sim.dat <- data.frame(sim.dat)

            sim.dat[,1] <- NA
            sim.dat[TIPS,1] <- phy$tip.label
            sim.dat[,2] <- c(data[,1], phy$node.label)
            sim.dat[,3] <- x
            
            if(tip.fog == "known"){
                for(i in TIPS){
                    sim.dat[i,3] <- rnorm(1,sim.dat[i,3], data[i,2])
                }
            }
        }else{
            sim.dat <- matrix(,ntips,3)
            sim.dat <- data.frame(sim.dat)

            sim.dat[,1] <- phy$tip.label
            sim.dat[,2] <- data[,1]
            sim.dat[,3] <- x[TIPS]
            
            if(tip.fog == "known"){
                for(i in TIPS){
                    sim.dat[i,3] <- rnorm(1,sim.dat[i,3], data[i,2])
                }
            }
        }

		colnames(sim.dat)<-c("Genus_species","Reg","X")
		if(tip.fog=="known"){
			sim.dat$tip.fog <- tip.fog_vector
		}
	}
	if(simmap.tree==TRUE){
		#This is annoying, but the second column has to be in there twice otherwise, error.
		if(tip.fog == "none"){
			data <- data.frame(data[,2], data[,2], row.names=data[,1])
		}
		if(tip.fog == "known"){
			data <- data.frame(data[,2], data[,3], row.names=data[,1])
			tip.fog_vector <- data[,2] #because we've shifted things over
		}
		if(is.numeric(tip.fog)){
			if(length(tip.fog) == length(phy$tip.label)){
				data <- data.frame(data[,2], tip.fog, row.names=data[,1])
				tip.fog_vector <- tip.fog
			}
			if(length(tip.fog)==1){
				data <- data.frame(data[,2], rep(tip.fog, length(phy$tip.label)), row.names=data[,1])
				tip.fog_vector <- rep(tip.fog, length(phy$tip.label))
			}
			tip.fog <- "known"
		}
		
		data <- data[phy$tip.label,]

		n=max(phy$edge[,1])
		ntips=length(phy$tip.label)

		k=length(colnames(phy$mapped.edge))

		regimeindex <- colnames(phy$mapped.edge)
		##Begins the construction of the edges matrix -- similar to the ouch format##
		#Makes a vector of absolute times in proportion of the total length of the tree
		branch.lengths=rep(0,(n-1))
		branch.lengths[(ntips+1):(n-1)]=branching.times(phy)[-1]/max(branching.times(phy))

		#Obtain root state and internal node labels
		root.edge.index <- which(phy$edge[,1] == ntips+1)
		root.state <- which(colnames(phy$mapped.edge)==names(phy$maps[[root.edge.index[2]]][1]))
		
		#New tree matrix to be used for subsetting regimes
		edges=cbind(c(1:(n-1)),phy$edge,MakeAgeTable(phy, root.age=root.age))
		if(scaleHeight==TRUE){
			edges[,4:5]<-edges[,4:5]/max(MakeAgeTable(phy, root.age=root.age))
			root.age <- max(MakeAgeTable(phy, root.age=root.age))
			phy$maps <- lapply(phy$maps, function(x) x/root.age)
            root.age = 1
		}
		edges=edges[sort.list(edges[,3]),]

		#Resort the edge matrix so that it looks like the original matrix order
		edges=edges[sort.list(edges[,1]),]

		oldregime=root.state
		oldtime=0

		alpha=alpha
		sigma=sqrt(sigma.sq)
		theta=theta

		n.cov=matrix(rep(0,n*n), n, n)
		nodecode=matrix(c(ntips+1,1),1,2)

		x <- matrix(0, n, 1)
		TIPS <- 1:ntips
		ROOT <- ntips + 1L
		x[ROOT,] <- theta0

		for(i in 1:length(edges[,1])){
			currentmap<-phy$maps[[i]]
			oldtime=edges[i,4]

			if(length(phy$maps[[i]])==1){
				regimeduration<-currentmap[1]
				newtime<-oldtime+regimeduration
				regimenumber<-which(colnames(phy$mapped.edge)==names(currentmap)[1])
				x[edges[i,3],]=x[edges[i,2],]*exp(-alpha[regimenumber]*(newtime-oldtime))+(theta[regimenumber])*(1-exp(-alpha[regimenumber]*(newtime-oldtime)))+sigma[regimenumber]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[regimenumber]*(newtime-oldtime)))/(2*alpha[regimenumber]))
			}
			if(length(phy$maps[[i]])>1){
				regimeduration<-currentmap[1]
				newtime<-oldtime+regimeduration
				regimenumber<-which(colnames(phy$mapped.edge)==names(currentmap)[1])
				x[edges[i,3],]=x[edges[i,2],]*exp(-alpha[regimenumber]*(newtime-oldtime))+(theta[regimenumber])*(1-exp(-alpha[regimenumber]*(newtime-oldtime)))+sigma[regimenumber]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[regimenumber]*(newtime-oldtime)))/(2*alpha[regimenumber]))
				oldtime<-newtime
				for (regimeindex in 2:length(currentmap)){
					regimeduration<-currentmap[regimeindex]
					newtime<-oldtime+regimeduration
					regimenumber<-which(colnames(phy$mapped.edge)==names(currentmap)[regimeindex])
					x[edges[i,3],]=x[edges[i,3],]*exp(-alpha[regimenumber]*(newtime-oldtime))+(theta[regimenumber])*(1-exp(-alpha[regimenumber]*(newtime-oldtime)))+sigma[regimenumber]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[regimenumber]*(newtime-oldtime)))/(2*alpha[regimenumber]))
                    oldtime<-newtime
                    newregime<-regimenumber
                }
            }
        }
        
        if(get.all == TRUE) {
            sim.dat <- matrix(,length(x),3)
            sim.dat <- data.frame(sim.dat)
            
            sim.dat[,1] <- NA
            sim.dat[TIPS,1] <- phy$tip.label
            sim.dat[,2] <- x
            
            if(tip.fog == "known"){
                for(i in TIPS){
                    sim.dat[i,2] <- rnorm(1, sim.dat[i,2], data[i,2])
                }
            }
        }else{
            sim.dat <- matrix(,ntips,2)
            sim.dat <- data.frame(sim.dat)
            
            sim.dat[,1] <- phy$tip.label
            sim.dat[,2] <- x[TIPS,]
			if(tip.fog == "known"){
                for(i in TIPS){
					print(sim.dat[i,2])
					print(data[i,2])
                    sim.dat[i,2] <- rnorm(1, sim.dat[i,2], data[i,2])
                }
            }
        }
        colnames(sim.dat)<-c("Genus_species","X")
		if(tip.fog=="known"){
			sim.dat$tip.fog <- tip.fog_vector
		}
    }
    sim.dat
}

