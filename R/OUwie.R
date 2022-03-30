#OUwie Master Controller

#written by Jeremy M. Beaulieu

#Fits the Ornstein-Uhlenbeck model of continuous characters evolving under discrete selective
#regimes. The input is a tree of class "phylo" that has the regimes as internal node labels
#and a trait file. The trait file must be in the following order: Species names, Regime, and
#continuous trait. Different models can be specified -- Brownian motion (BM), multiple rate BM (BMS)
#global OU (OU1), multiple regime OU (OUM), multiple sigmas (OUMV), multiple alphas (OUMA),
#and the multiple alphas and sigmas (OUMVA).

OUwie <- function(phy, data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA", "TrendyM", "TrendyMS"), simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, get.root.theta=FALSE, shift.point=0.5, clade=NULL, mserr="none", starting.vals=NULL, check.identify=TRUE, algorithm=c("invert", "three.point"), diagn=FALSE, quiet=FALSE, warn=TRUE, lb = NULL, ub = NULL, opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5)){

    if(length(algorithm) == 2){
        algorithm = "invert"
        warning("An algorithm was not specified. Defaulting to computing the determinant and inversion of the vcv.", call.=FALSE, immediate.=TRUE)
    }
	
	if(algorithm == "fast") {
		algorithm = "three.point"	
	}
	
	if(algorithm == "slow") {
		algorithm = "invert"	
	}	

    if(model=="BMS" & root.station==TRUE){
        warning("By setting root.station=TRUE, you have specified the group means model of Thomas et al. 2006", call.=FALSE, immediate.=TRUE)
        get.root.theta = FALSE
    }

    if(is.factor(data[,3])==TRUE){
        stop("Check the format of the data column. It's reading as a factor.", call. = FALSE)
    }

    if(!is.null(starting.vals[1])){
        if(model == "OU1" | model == "OUM" | model == "OUMV" | model == "OUMA" | model == "OUMVA"){
            if(length(starting.vals)<2){
                stop("You only supplied one starting value. For OU models you need to supply two", call. = FALSE)
            }
        }
    }
    
    if(model == "OUMV" |  model == "OUMA" |  model == "OUMVA"){
        if(root.station == TRUE){
            stop("Assuming stationarity at the root is no longer allowed under these models. Try OU1 and OUM.", call. = FALSE)
        }
    }

    if(model == "BM1" |  model == "BMS" | model == "TrendyM" | model == "TrendyMS"){
        if(root.station == FALSE){
            get.root.theta = TRUE
        }
    }

    if(diagn == TRUE & algorithm == "three.point"){
        diagn=FALSE
        warning("Turning off internal diagnostic, not implemented for the three-point algorithm.", call. = FALSE, immediate.=TRUE)
    }
    
    if(check.identify == TRUE){
            identifiable <- check.identify(phy=phy, data=data, simmap.tree=simmap.tree, quiet=TRUE)
            if(identifiable == 0){
                warning("The supplied regime painting may be unidentifiable for the regime painting. All regimes form connected subtrees.", call. = FALSE, immediate.=TRUE)
            }
    }


    #Makes sure the data is in the same order as the tip labels
    if(mserr == "none"){
		data <- data.frame(data[,2], data[,3], row.names=data[,1])
		data <- data[phy$tip.label,]
		tip.states.cp <- factor(data[,1]) # fixes a cosmetic bug when internal states don't match external
	}
	if(mserr=="known"){
		if(!dim(data)[2]==4){
			stop("You specified measurement error should be incorporated, but this information is missing", call. = FALSE)
		}
		else{
            if(is.factor(data[,4]) == TRUE){
                stop("Check the format of the measurement error column. It's reading as a factor.", call. = FALSE)
            }else{
                data <- data.frame(data[,2], data[,3], data[,4], row.names=data[,1])
                data <- data[phy$tip.label,]
                tip.states.cp <- factor(data[,1])
            }
		}
	}
  
  # if(algorithm == "three.point" & scaleHeight == FALSE){
  #   warning("It is recommended that you set scaleHeight to TRUE if using the three-point algorithm.", call. = FALSE, immediate.=TRUE)
  # }

    #Values to be used throughout
	n <- max(phy$edge[,1])
	ntips <- length(phy$tip.label)
	#Will label the clade of interest if the user so chooses:
	if(!is.null(clade) & simmap.tree==FALSE){
		node <- mrca(phy)[clade[1],clade[2]]
		int <- c(node,Descendants(phy,node, "all"))
		tips <- int[int<ntips]
		data[,1] <- 1
		data[tips,1] <- 2
		int <- int[int>ntips]
		phy$node.label <- rep(1,phy$Nnode)
		pp <- int-length(phy$tip.label)
		phy$node.label[pp] <- 2
	}
	if (is.character(model)) {
		if (model == "BM1"| model == "OU1"){
			simmap.tree=FALSE
		}
	}
	if(simmap.tree==TRUE){
		k <- length(colnames(phy$mapped.edge))
		tot.states <- factor(colnames(phy$mapped.edge))
		tip.states <- factor(data[,1])
		data[,1] <- as.numeric(tip.states)

		#Obtains the state at the root
		root.edge.index <- which(phy$edge[,1] == ntips+1)
		root.state <- which(colnames(phy$mapped.edge)==names(phy$maps[[root.edge.index[2]]][1]))
		##Begins the construction of the edges matrix -- similar to the ouch format##
		edges=cbind(c(1:(n-1)),phy$edge, MakeAgeTable(phy, root.age=root.age))

		if(scaleHeight==TRUE){
            Tmax <- max(MakeAgeTable(phy, root.age=root.age))
			edges[,4:5]<-edges[,4:5]/Tmax
            root.age = 1
		}
		edges=edges[sort.list(edges[,3]),]

		#Resort the edge matrix so that it looks like the original matrix order
		edges=edges[sort.list(edges[,1]),]
	}
	if(simmap.tree==FALSE){

		#A boolean for whether the root theta should be estimated -- default is that it should be.
		if (is.character(model)) {
			if (model == "BM1"| model == "OU1"){
				##Begins the construction of the edges matrix -- similar to the ouch format##
				#Makes a vector of absolute times in proportion of the total length of the tree
                k <- 1
                phy$node.label <- rep(1, Nnode(phy))
                tip.states <- factor(rep(1, length(data[,1])))
                int.states <- factor(phy$node.label)
                tot.states <- factor(c(phy$node.label,tip.states))

                #Obtain root state -- for both models assume the root state to be 1 since no other state is used even if provided in the tree
				root.state <- 1
				#New tree matrix to be used for subsetting regimes
				edges=cbind(c(1:(n-1)),phy$edge,MakeAgeTable(phy, root.age=root.age))
				if(scaleHeight==TRUE){
                    Tmax <- max(MakeAgeTable(phy, root.age=root.age))
                    edges[,4:5]<-edges[,4:5]/Tmax
                    root.age <- 1
				}
				edges=edges[sort.list(edges[,3]),]

				regime <- matrix(0,nrow=length(edges[,1]),ncol=k)
				regime[,1] <- 1
				if (k>1) {
					regime[,2:k] <- 0
				}
				edges=cbind(edges,regime)
			}
			else{
                #Obtain a a list of all the regime states. This is a solution for instances when tip states and
                #the internal nodes are not of equal length:
                tot.states <- factor(c(phy$node.label,as.character(data[,1])))
                k <- length(levels(tot.states))
                int.states <- factor(phy$node.label)
                phy$node.label <- as.numeric(int.states)
                tip.states <- factor(data[,1])
                data[,1] <- as.numeric(tip.states)
				#Obtain root state and internal node labels
				root.state <- phy$node.label[1]
				int.state <- phy$node.label[-1]
				#New tree matrix to be used for subsetting regimes
                edges <- cbind(c(1:(n-1)),phy$edge,MakeAgeTable(phy, root.age=root.age))
				if(scaleHeight==TRUE){
                    Tmax <- max(MakeAgeTable(phy, root.age=root.age))
                    edges[,4:5]<-edges[,4:5]/Tmax
                    root.age <- 1
				}
				edges=edges[sort.list(edges[,3]),]

				mm <- c(data[,1],int.state)
				regime <- matrix(0,nrow=length(mm),ncol=k)
				#Generates an indicator matrix from the regime vector
				for (i in 1:length(mm)) {
					regime[i,mm[i]] <- 1
				}
				#Finishes the edges matrix
				edges=cbind(edges,regime)
			}
		}
		#Resort the edge matrix so that it looks like the original matrix order
		edges <- edges[sort.list(edges[,1]),]
	}

    #Initializes an integer specifying whether or not the trendy model is to be invoked. Assume 0, unless otherwise stated later:
	trendy <- 0
	
    #Data:
    if(algorithm == "three.point"){
        x <- data[,2]
        names(x) <- rownames(data)
        tip.paths <- lapply(1:length(data[,2]), function(x) getPathToRoot(phy, x))
    }else{
        x <- as.matrix(data[,2])
    }

    #Matches the model with the appropriate parameter matrix structure
	if (is.character(model)) {
		index.mat <- matrix(0, 2,k)
        if(algorithm == "three.point"){
            Rate.mat <- matrix(1, 3, k)
        }else{
            Rate.mat <- matrix(1, 2, k)
        }
		if (model == "BM1"){
			np <- 1
            index.mat[1,1:k] <- NA
			index.mat[2,1:k] <- 1
            if(algorithm == "three.point"){
                index.mat <- rbind(index.mat, rep(NA,k))
            }
			param.count <- np+1
		}
		#The group mean model of Thomas et al, is trivial: set root.station to be TRUE:
		if (model == "BMS"){
			np <- k
            index.mat[1,1:k] <- NA
			index.mat[2,1:k] <- 1:np
			if(root.station==TRUE){
                if(algorithm == "three.point"){
                    max.par.so.far <- max(index.mat)
                    index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
                }
				param.count <- np+k
			}
            if(root.station==FALSE){
                if(algorithm == "three.point"){
                    index.mat <- rbind(index.mat, rep(NA,k))
                }
                param.count <- np+1
            }
		}
		if (model == "OU1"){
			np <- 2
			index.mat[1,1:k] <- 1
			index.mat[2,1:k] <- 2
            if(algorithm == "three.point"){
                max.par.so.far <- max(index.mat)
                index.mat <- rbind(index.mat, rep(3,k))
            }
            param.count <- np + 1
            if(get.root.theta == TRUE){
                param.count <- param.count + 1
            }
		}
		if (model == "OUM"){
            np <- 2
			index.mat[1,1:k] <- 1
			index.mat[2,1:k] <- 2
            if(algorithm == "three.point"){
                max.par.so.far <- max(index.mat)
                index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
            }
            param.count <- np + k
            if(get.root.theta == TRUE){
                param.count <- param.count + 1
            }
		}
		if (model == "OUMV") {
			np <- k+1
			index.mat[1,1:k] <- 1
			index.mat[2,1:k] <- 2:(k+1)
            if(algorithm == "three.point"){
                max.par.so.far <- max(index.mat)
                index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
            }
            param.count <- np + k
            if(get.root.theta == TRUE){
                param.count <- param.count + 1
            }
		}
		if (model == "OUMA") {
			np <- k+1
			index.mat[1,1:k] <- 1:k
			index.mat[2,1:k] <- k+1
            if(algorithm == "three.point"){
                max.par.so.far <- max(index.mat)
                index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
            }
            param.count <- np+k
            if(get.root.theta == TRUE){
                param.count <- param.count + 1
            }
		}
		if (model == "OUMVA") {
			np <- k*2
            index <- matrix(TRUE,2,k)
			index.mat[index] <- 1:(k*2)
            if(algorithm == "three.point"){
                max.par.so.far <- max(index.mat)
                index.mat <- rbind(index.mat, (max.par.so.far + 1):(max.par.so.far + k))
            }
            param.count <- np+k
            if(get.root.theta == TRUE){
                param.count <- param.count + 1
            }
		}
		if (model == "TrendyM"){
			index.mat <- matrix(0,3,k)
			np <- k+2 #We have k regimes, but 1 sigma, plus the root value
            index.mat[1,1:k] <- np+1
			index.mat[2,1:k] <- 1
			index.mat[3,1:k] <- 2:(np-1)
			param.count<-np
			root.station=TRUE
			trendy=1
		}
		if (model == "TrendyMS"){
			index.mat <- matrix(0,3,k)
			np <- (k*2)+1 #We have k regimes and assume k trends and k sigmas, plus the root value
            index.mat[1,1:k] <- np+1
			index.mat[2,1:k] <- 1:k
			index.mat[3,1:k] <- (k+1):(np-1)
			param.count <- np
			root.station=FALSE
			trendy=1
		}
	}

	if(model == "BM1" | model == "BMS"){
	  if(algorithm == "three.point"){
	    index.mat[is.na(index.mat)] <- param.count + 1
	  }else{
	    index.mat[is.na(index.mat)] <- param.count
	  }
	}else{
	  index.mat[is.na(index.mat)] <- param.count + 1
	}	
	
	if(warn==TRUE){
		if(param.count > (ntips/10)){
			warning("You might not have enough data to fit this model well", call.=FALSE, immediate.=TRUE)
		}
	}

	#Likelihood function for estimating model parameters
	dev <- function(p, index.mat, edges, mserr, trendy, get.root.theta){
		if(trendy == 0){
            
            p <- exp(p)
            Rate.mat[] <- c(p, 1e-10)[index.mat]
            if(get.root.theta == TRUE){
                W <- weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=FALSE, shift.point=shift.point)
            }else{
                W <- weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=TRUE, shift.point=shift.point)
            }
            
            if(algorithm == "invert"){
                N <- length(x[,1])
                V <- varcov.ou(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=root.station, shift.point=shift.point)
                
                if (any(is.nan(diag(V))) || any(is.infinite(diag(V)))) return(1000000)
                if(mserr=="known"){
                    diag(V) <- diag(V)+(data[,3]^2)
                }
                
                theta <- Inf
                try(theta <- pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%x, silent=TRUE)
                if(any(theta==Inf)){
                    return(10000000)
                }
                DET <- sum(log(abs(Re(diag(qr(V)$qr)))))
                #However, sometimes this fails (not sure yet why) so I just toggle between this and another approach:
                if(!is.finite(DET)){
                    DET <- determinant(V, logarithm=TRUE)
                    logl <- -.5*(t(x-W%*%theta)%*%pseudoinverse(V)%*%(x-W%*%theta))-.5*as.numeric(DET$modulus)-.5*(N*log(2*pi))
                }else{
                    logl <- -.5*(t(x-W%*%theta)%*%pseudoinverse(V)%*%(x-W%*%theta))-.5*as.numeric(DET)-.5*(N*log(2*pi))
                }
            }
            if(algorithm == "three.point"){
                pars <- matrix(c(Rate.mat[3,], Rate.mat[2,], Rate.mat[1,]), dim(Rate.mat)[2], 3, 
                               dimnames = list(levels(factor(as.numeric(tot.states))), c("opt", "sig", "alp")))
                if(get.root.theta == TRUE){
                    root.par.index <- length(p)
                    theta0 <- p[root.par.index]
                    expected.vals <- colSums(t(W) * c(theta0, pars[,1]))
                    names(expected.vals) <- phy$tip.label
                }else{
                    expected.vals <- colSums(t(W) * pars[,1])
                    names(expected.vals) <- phy$tip.label
                }
                transformed.tree <- transformPhy(phy, map, pars, tip.paths)
                if(mserr=="known"){
                  TIPS <- transformed.tree$tree$edge[,2] <= length(transformed.tree$tree$tip.label)
                  transformed.tree$tree$edge.length[TIPS] <- transformed.tree$tree$edge.length[TIPS] + (data[,3]^2/transformed.tree$diag/transformed.tree$diag)
                }
                comp <- NA
                try(comp <- phylolm::three.point.compute(transformed.tree$tree, x, expected.vals, transformed.tree$diag), silent=TRUE)
                if(is.na(comp[1])){
                    return(10000000)
                }else{
                    logl <- -as.numeric(Ntip(phy) * log(2 * pi) + comp$logd + (comp$PP - 2 * comp$QP + comp$QQ))/2
                }
            }
		}else{
            p <- exp(p)
            Rate.mat[] <- c(p, 1e-10)[index.mat]
            N <- length(x[,1])
            root.par.index <- length(p)
            V <- varcov.ou(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=root.station, shift.point=shift.point)
            if (any(is.nan(diag(V))) || any(is.infinite(diag(V)))) return(1000000)
            if(mserr=="known"){
                diag(V)<-diag(V)+(data[,3]^2)
            }
			E_a <- Inf
			try(E_a <- expected.trendy(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, root.value=p[root.par.index], shift.point=shift.point))
			if(any(E_a==Inf)){
				return(10000000)
			}
			DET <- sum(log(abs(Re(diag(qr(V)$qr)))))
			#However, sometimes this fails (not sure yet why) so I just toggle between this and another approach:
			if(!is.finite(DET)){
				DET <- determinant(V, logarithm=TRUE)
				logl <- -.5*(t(x-E_a)%*%pseudoinverse(V)%*%(x-E_a))-.5*as.numeric(DET$modulus)-.5*(N*log(2*pi))
			}else{
				logl <- -.5*(t(x-E_a)%*%pseudoinverse(V)%*%(x-E_a))-.5*as.numeric(DET)-.5*(N*log(2*pi))
			}
		}
        
		#When the model includes alpha, the values of V can get too small, the modulus does not seem correct and the loglik becomes unstable. This is one solution:
		if(!is.finite(logl)){
			return(10000000)
		}
		return(-logl)
	}

	if(quiet==FALSE){
		cat("Initializing...", "\n")
	}
  
	# upper and lower bounds
	if(algorithm == "three.point"){
	  if(max(index.mat[1,]) > max(index.mat[2,])){
	    k.alpha <- 0
	    k.sigma <- length(unique(index.mat[2,]))
	  }else{
	    k.alpha <- length(unique(index.mat[1,]))
	    k.sigma <- length(unique(index.mat[2,]))
	  }
	  if(is.null(lb)){
	    lb.alpha <- 1e-9
	    lb.sigma <- 1e-9
	    lower = c(rep(log(lb.alpha), k.alpha), rep(log(lb.sigma), k.sigma))
	    lb <- 1e-9
	  }else{
	    lower = c(rep(log(lb[1]), k.alpha), rep(log(lb[2]), k.sigma))
	    ub <- ub[3] # theta's are added later
	  }
	  if(is.null(ub)){
	    ub.alpha <- 100
	    ub.sigma <- 100
	    upper = c(rep(log(ub.alpha), k.alpha), rep(log(ub.sigma), k.sigma))
	    ub <- 100
	  }else{
	    upper = c(rep(log(ub[1]), k.alpha), rep(log(ub[2]), k.sigma))
	    ub <- ub[3] # theta's are added later
	  }
	}else{
	  if(is.null(lb)){
	    lb <- 1e-9
	  }
	  if(is.null(ub)){
	    ub <- 100
	  }
	  lower = rep(log(lb), np)
	  upper = rep(log(ub), np)	
	}

	# for any downstream usage, we default to the old treatment of lb and ub
	
    if(algorithm == "three.point"){
        if(simmap.tree == FALSE){
            map <- getMapFromNode(phy, tip.states, int.states, shift.point)
            if(scaleHeight==TRUE){
              map <- lapply(map, function(x) x/Tmax)
            }
        }else{
            map <- phy$maps
            if(scaleHeight==TRUE){
              map <- lapply(map, function(x) x/Tmax)
              phy$maps <- map
            }
        }
    }

	if(scaleHeight==TRUE){
	  phy$edge.length <- phy$edge.length/Tmax
	  Tmax <- 1
	}else{
	  Tmax <- max(branching.times(phy))
	}
	
	if(model == "OU1" | model == "OUM" | model == "OUMV" | model == "OUMA" | model == "OUMVA"){
        #Initial value for alpha is just the half life based on the entire length of the tree:
        if(is.null(starting.vals)){
            init.ip <- c(log(2)/Tmax, mean(pic(x, phy)^2))
        }else{
            init.ip <- starting.vals
        }
		if(model=="OU1"){
			ip <- init.ip
		}
		if(model=="OUMV" | model=="OUMA" | model=="OUM"){
			ip <- c(rep(init.ip[1],length(unique(index.mat[1,]))), rep(init.ip[2],length(unique(index.mat[2,]))))
		}
		if(model=="OUMVA"){
			ip <- c(rep(init.ip,k))
		}
        if(algorithm == "three.point"){
            if(model == "OU1"){
                ip <- c(ip, mean(x))
            }else{
                means.by.regime <- with(data, tapply(data[,2], data[,1], mean))
                if(length(means.by.regime) < length(levels(tot.states))){
                  means.by.regime <- rep(mean(data[,2]), length(levels(tot.states)))
                }
                names(means.by.regime) <- NULL
                ip <- c(ip, means.by.regime)
            }
            lower <- c(lower, rep(log(lb), k))
            upper <- c(upper, rep(log(ub), k))
            if(get.root.theta == TRUE){
                ip <- c(ip, means.by.regime[root.state])
                lower <- c(lower, log(lb))
                upper <- c(upper, log(ub))
            }
        }
		if(quiet==FALSE){
			cat("Finished. Begin thorough search...", "\n")
		}
		out <- nloptr(x0=log(ip), eval_f=dev, lb=lower, ub=upper, opts=opts, index.mat=index.mat, edges=edges, mserr=mserr, trendy=trendy, get.root.theta=get.root.theta)
	}
	else{
        if(is.null(starting.vals)){
            sig <- mean(pic(x, phy)^2)
        }else{
            sig <- starting.vals
        }
		#####################
	  if(model=="BMS"){
	    ip <- rep(sig, k)
	    if(get.root.theta == TRUE){
	      if(algorithm == "three.point"){
	        ip <- c(ip, mean(x))
	        lower <- c(lower, log(lb))
	        upper <- c(upper, log(ub))
	      }
	    }
	  }
	  else{
	    if(model=="BM1"){
	      ip <- sig
	      if(get.root.theta == TRUE){
	        if(algorithm == "three.point"){
	          ip <- c(ip, mean(x))
	          lower <- c(lower, log(lb))
	          upper <- c(upper, log(ub))
	        }
	      }
	    }
	    if(model=="TrendyM"){
				#We assume that the starting trend values are zero:
				ip <- c(sig, rep(exp(-20), param.count-2), mean(x))
				lower <- c(lb,rep(log(lb), param.count-2), log(lb))
				upper <- c(ub,rep(log(ub), param.count-2), log(ub))
			}
			if(model == "TrendyMS"){
				#We assume that the starting trend values are zero:
				ip <- c(rep(sig, k), rep(exp(-20), (np-1) - k), mean(x))
				lower <- c(rep(lb,k), rep(log(lb), (np-1) - k), log(lb))
				upper <- c(rep(ub,k), rep(log(ub), (np-1) - k), log(ub))
			}
		}
		if(quiet==FALSE){
			cat("Finished. Begin thorough search...", "\n")
		}
		out = nloptr(x0=log(ip), eval_f=dev, lb=fix_lower(lower, log(ip)), ub=fix_upper(upper, log(ip)), opts=opts, index.mat=index.mat, edges=edges, mserr=mserr, trendy=trendy, get.root.theta=get.root.theta)
	}
	
    loglik <- -out$objective
	out$solution <- exp(out$solution)
	out$new.starting <- out$solution # note: back in regular param space, not log space used in the search
	
    #Takes estimated parameters from dev and calculates theta for each regime:
	dev.theta <- function(p, index.mat, edges=edges, mserr=mserr){
		Rate.mat[] <- c(p, 1e-10)[index.mat]
		N <- length(x[,1])
		V <- varcov.ou(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=root.station, shift.point=shift.point)
        
        if(get.root.theta == TRUE){
            W <- weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=FALSE, shift.point=shift.point)
        }else{
            W <- weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=TRUE, shift.point=shift.point)
        }
        
        if(mserr=="known"){
            diag(V) <- diag(V)+(data[,3]^2)
		}
		theta <- pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%x
        
		#Standard error of theta -- uses pseudoinverse to overcome singularity issues
		se <- sqrt(diag(pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)))
		#Joins the vector of thetas with the vector of standard errors into a 2 column matrix for easy extraction at the summary stage
		theta.est <- cbind(theta,se)
        return(theta.est)
	}

	#Informs the user that the summarization has begun, output model-dependent and dependent on whether the root theta is to be estimated
	if(quiet==FALSE){
		cat("Finished. Summarizing results.", "\n")
	}

	if(model == "TrendyM" | model == "TrendyMS"){
		theta<-NULL
		root.est <- out$solution[length(out$solution)]
		theta <- matrix(root.est, k+1, 2)
		theta[,2] <- 0
	}else{
        if(algorithm == "invert"){
            theta <- dev.theta(out$solution, index.mat, edges, mserr)
        }else{
            if(model == "BM1" | model == "BMS"){
                max.pars <- max(index.mat)
                index.mat[3,] <- max.pars - 1
            }
        }
	}

	#Calculates the Hessian for use in calculating standard errors and whether the maximum likelihood solution was found
	if(diagn==TRUE){
	  data[,1] <- tip.states.cp # fixes a cosmetic bug in cases where internal states don't match tip states
		h <- hessian(x=log(out$solution), func=dev, index.mat=index.mat, edges=edges, mserr=mserr, trendy=trendy, get.root.theta=get.root.theta)
		#Using the corpcor package here to overcome possible NAs with calculating the SE
		solution <- matrix(out$solution[index.mat], dim(index.mat))
		solution.se <- matrix(sqrt(diag(pseudoinverse(h)))[index.mat], dim(index.mat))
        if(algorithm == "invert"){
            rownames(solution) <- rownames(solution.se) <- rownames(index.mat) <- c("alpha","sigma.sq")
        }
        if(algorithm == "three.point"){
            rownames(solution) <- rownames(solution.se) <- rownames(index.mat) <- c("alpha","sigma.sq", "theta")
        }
		if(simmap.tree==FALSE){
			colnames(solution) <- colnames(solution.se) <- levels(tot.states)
		}
		if(simmap.tree==TRUE){
			colnames(solution) <- colnames(solution.se) <- c(colnames(phy$mapped.edge))
		}
		#Eigendecomposition of the Hessian to assess reliability of likelihood estimates
		hess.eig<-eigen(h,symmetric=TRUE)
		#If eigenvect is TRUE then the eigenvector and index matrix will appear in the list of objects
		eigval<-signif(hess.eig$values,2)
		eigvect<-round(hess.eig$vectors, 2)
        if(get.root.theta == TRUE){
            if(model == "BM1" | model == "BMS"){
                W <- matrix(0, Ntip(phy), k + 1)
                W[,1] <- 1
            }else{
                W <- weight.mat(phy, edges, Rate.mat=solution, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=FALSE, shift.point=shift.point)
            }
            rates <- cbind(NA, solution, NA)
            if(algorithm == "invert"){
                thetas <- cbind(t(theta[,1]), NA)
                rates <- rbind(rates, thetas)
            }else{
                theta.est <- solution[3,]
                se <- rep(NA,length(theta))
                theta <- cbind(theta.est,se)
            }
            weights <- cbind(W, data[,1])
            regime.weights <- rbind(rates, weights)
            rownames(regime.weights) <- c("alpha", "sigma.sq", "theta", phy$tip.label)
            colnames(regime.weights) <- c("Root", levels(tot.states), "Tip_regime")
        }else{
            if(model=="TrendyM" | model=="TrendyMS"){
                regime.weights = NULL
            }else{
                W <- weight.mat(phy, edges, Rate.mat=solution, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=TRUE, shift.point=shift.point)
                rates <- cbind(solution, NA)
                if(algorithm == "invert"){
                    thetas <- cbind(t(theta[,1]), NA)
                    rates <- rbind(rates, thetas)
                }else{
                    theta.est <- solution[3,]
                    se <- rep(NA,length(theta))
                    theta <- cbind(theta.est,se)
                }
                weights <- cbind(W, data[,1])
                regime.weights <- rbind(rates, weights)
                rownames(regime.weights) <- c("alpha", "sigma.sq", "theta", phy$tip.label)
                colnames(regime.weights) <- c(levels(tot.states), "Tip_regime")
            }
        }

		obj = list(loglik = loglik, AIC = -2*loglik+2*param.count, AICc=-2*loglik+(2*param.count*(ntips/(ntips-param.count-1))), BIC = -2*loglik + log(ntips) * param.count, model=model, param.count=param.count, solution=solution, theta=theta, solution.se=solution.se, tot.states=tot.states, index.mat=index.mat, simmap.tree=simmap.tree, root.age=root.age, shift.point=shift.point, opts=opts, data=data, phy=phy, root.station=root.station, scaleHeight=scaleHeight, starting.vals=starting.vals, lb=lower, ub=upper, iterations=out$iterations, get.root.theta=get.root.theta, regime.weights=regime.weights, eigval=eigval, eigvect=eigvect, new.start=out$new.start, algorithm=algorithm)

	}
    
	if(diagn==FALSE){
	  data[,1] <- tip.states.cp # fixes a cosmetic bug in cases where internal states don't match tip states
		solution <- matrix(out$solution[index.mat], dim(index.mat))
		if(model=="TrendyM" | model=="TrendyMS"){
			rownames(solution) <- rownames(index.mat) <- c("alpha","sigma.sq","trend")
		}else{
            if(algorithm == "three.point"){
                rownames(solution) <- rownames(index.mat) <- c("alpha","sigma.sq", "theta")
            }else{
                rownames(solution) <- rownames(index.mat) <- c("alpha","sigma.sq")
            }
		}
		if(simmap.tree==FALSE){
			colnames(solution) <- levels(tot.states)
		}
		if(simmap.tree==TRUE){
			colnames(solution) <- c(colnames(phy$mapped.edge))
		}
		if(mserr=="est"){
			mserr.est<-out$solution[length(out$solution)]
			param.count<-param.count+1
		}
		else{
			mserr.est<-NULL
		}
        
        if(get.root.theta == TRUE){
            if(model == "BM1" | model == "BMS"){
                W <- matrix(0, Ntip(phy), k+1)
                W[,1] <- 1
            }else{
                W <- weight.mat(phy, edges, Rate.mat=solution, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=FALSE, shift.point=shift.point)
            }
            rates <- cbind(NA, solution, NA)
            if(algorithm == "invert"){
                thetas <- cbind(t(theta[,1]), NA)
                rates <- rbind(rates, thetas)
            }else{
                theta.est <- solution[3,]
                se <- rep(NA,length(theta.est))
                theta <- cbind(theta.est,se)
            }
            weights <- cbind(W, data[,1])
            regime.weights <- rbind(rates, weights)
            rownames(regime.weights) <- c("alpha", "sigma.sq", "theta", phy$tip.label)
            colnames(regime.weights) <- c("Root", levels(tot.states), "Tip_regime")
        }else{
            if(model=="TrendyM" | model=="TrendyMS"){
                regime.weights = NULL
            }else{
                W <- weight.mat(phy, edges, Rate.mat=solution, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=TRUE, shift.point=shift.point)
                rates <- cbind(solution, NA)
                if(algorithm == "invert"){
                    thetas <- cbind(t(theta[,1]), NA)
                    rates <- rbind(rates, thetas)
                }else{
                    theta.est <- solution[3,]
                    se <- rep(NA,length(theta.est))
                    theta <- cbind(theta.est,se)
                }
                weights <- cbind(W, data[,1])
                regime.weights <- rbind(rates, weights)
                rownames(regime.weights) <- c("alpha", "sigma.sq", "theta", phy$tip.label)
                colnames(regime.weights) <- c(levels(tot.states), "Tip_regime")
            }
        }

        obj = list(loglik = loglik, AIC = -2*loglik+2*param.count,AICc=-2*loglik+(2*param.count*(ntips/(ntips-param.count-1))), BIC=-2*loglik + log(ntips) * param.count, model=model, param.count=param.count, solution=solution, theta=theta, tot.states=tot.states, index.mat=index.mat, simmap.tree=simmap.tree, root.age=root.age, shift.point=shift.point, opts=opts, data=data, phy=phy, root.station=root.station, scaleHeight=scaleHeight, starting.vals=starting.vals, lb=lower, ub=upper, iterations=out$iterations, get.root.theta=get.root.theta, regime.weights=regime.weights, new.start=out$new.start, algorithm=algorithm)
	}
	class(obj)<-"OUwie"
	return(obj)
}


print.OUwie<-function(x, ...){

	ntips <- Ntip(x$phy)
	output <- data.frame(x$loglik,x$AIC,x$AICc,x$BIC,x$model,ntips, row.names="")
	names(output) <- c("lnL","AIC","AICc","BIC","model","ntax")
	cat("\nFit\n")
	print(output)
	cat("\n")

	if (is.character(x$model)) {
		if (x$model == "BM1" | x$model == "BMS" | x$model == "TrendyM" | x$model == "TrendyMS"){
            param.est<- x$solution[1:2,]
            if(x$root.station==FALSE){
                theta.mat <- matrix(t(x$theta[1,]), 2, length(levels(x$tot.states)))
            }
			else{
				theta.mat<-matrix(t(x$theta), 2, length(levels(x$tot.states)))
			}
			rownames(theta.mat)<-c("estimate", "se")
			if(x$simmap.tree==FALSE){
				colnames(theta.mat) <- levels(x$tot.states)
			}
			if(x$simmap.tree==TRUE){
				colnames(theta.mat) <- c(colnames(x$phy$mapped.edge))
			}
			cat("Rates\n")
            print(param.est)
			cat("\n")
			cat("Optima\n")
			print(theta.mat)
			cat("\n")
		}
        if (x$get.root.theta == FALSE){
            if (x$model == "OU1" | x$model == "OUM"| x$model == "OUMV"| x$model == "OUMA" | x$model == "OUMVA"){
                param.est<- x$solution[1:2,]
                theta.mat<-matrix(t(x$theta), 2, length(levels(x$tot.states)))
                rownames(theta.mat)<-c("estimate", "se")
                if(x$simmap.tree==FALSE){
                    colnames(theta.mat)<- levels(x$tot.states)
                }
                if(x$simmap.tree==TRUE){
                    colnames(theta.mat) <- c(colnames(x$phy$mapped.edge))
                }
                cat("\nRates\n")
                print(param.est)
                cat("\n")
                cat("Optima\n")
                print(theta.mat)
                cat("\n")
                cat("\nHalf life (another way of reporting alpha)\n")
                if(x$model == "OU1"){
                    print(log(2)/param.est[1])
                }else{
                    print(log(2)/param.est['alpha',])
                }
                cat("\n")
            }
        }
        if (x$get.root.theta == TRUE){
            if (x$model == "OU1" | x$model == "OUM"| x$model == "OUMV"| x$model == "OUMA" | x$model == "OUMVA"){
                param.est<- x$solution[1:2,]
                theta.mat<-matrix(t(x$theta), 2, length(levels(x$tot.states))+1)
                rownames(theta.mat)<-c("estimate", "se")
                if(x$simmap.tree==FALSE){
                    colnames(theta.mat)<-c("root", levels(x$tot.states))
                }
                if(x$simmap.tree==TRUE){
                    colnames(theta.mat)<-c("root", colnames(x$phy$mapped.edge))
                }
                cat("\nRates\n")
                print(param.est)
                cat("\n")
                cat("Optima\n")
                print(theta.mat)
                cat("\n")
                cat("\nHalf life (another way of reporting alpha)\n")
                if(x$model == "OU1"){
                    print(log(2)/param.est[1])
                }else{
                    print(log(2)/param.est['alpha',])
                }
                cat("\n")
            }
        }
    }
	if(any(x$eigval<0)){
		index.matrix <- x$index.mat
		if(x$simmap.tree==FALSE){
			colnames(index.matrix) <- levels(x$tot.states)
		}
		if(x$simmap.tree==TRUE){
			colnames(index.matrix) <- colnames(x$phy$mapped.edge)
		}
		#If any eigenvalue is less than 0 then the solution is not the maximum likelihood solution
		if (any(x$eigval<0)) {
			cat("The objective function may be at a saddle point -- check eigenvectors or try a simpler model", "\n")
		}
	} else if (abs(x$loglik)==10000000) {
		cat("Error in calculating likeliood so an inaccurate placeholder value was returned: search likely never went beyond its starting value and the likelihood and AICc values are incorrect", "\n")
	} else {
		cat("Arrived at a reliable solution","\n")
	}
}



MakeAgeTable <- function(phy, root.age=NULL){
    if(is.null(root.age)){
        node.ages <- dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))
    }else{
        node.ages <- dateNodes(phy, rootAge=root.age)
    }
    max.age <- max(node.ages)
    table.ages <- matrix(0, dim(phy$edge)[1], 2)
    for(row.index in 1:dim(phy$edge)[1]){
        table.ages[row.index,1] <- max.age - node.ages[phy$edge[row.index,1]]
        table.ages[row.index,2] <- max.age - node.ages[phy$edge[row.index,2]]
    }
    return(table.ages)
}
