#Simple bootstrap function for assessing confidence in OUwie estimates

#written by Jeremy M. Beaulieu

OUwie.boot <- function(phy, data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"), nboot=100, alpha, sigma.sq, theta, theta0, simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, get.root.theta=FALSE, shift.point=0.5, clade=NULL, tip.fog="none", algorithm=c("invert", "three.point"), diagn=FALSE, quiet=TRUE, warn=FALSE){
	
    if(length(algorithm) == 2){
        algorithm = "invert"
        warning("An algorithm was not specified. Defaulting to computing the determinant and inversion of the vcv.", call.=FALSE, immediate.=TRUE)
    }

    if(is.null(root.age)){
        if(any(branching.times(phy)<0)){
            stop("Looks like your tree is producing negative branching times. Must input known root age of tree.", .call=FALSE)
        }
    }

	#the matrix we are building to store the results:
	res <- c()
	#if alpha is NA set to a really small number -- this is only relevant for BM models:
	alpha[is.na(alpha)] <- 1e-10

    cat("Beginning parametric bootstrap -- performing", nboot, "replicates", "\n")
	
	for(i in 1:nboot){
		tmp.phy<-phy
		#if(tip.fog=="known"){
			#else{
				#Now lengthen the terminal branches to reflect the intraspecific variation at the tips:
				#terminals <- tmp.phy$edge[,2] <= Ntip(tmp.phy)
				#terminal.edges <- tmp.phy$edge.length[terminals]
				#tmp.phy$edge.length[terminals] <- terminal.edges + data[,4]
			#}
		#}
		#This calls the OUwie simulator and simulates datasets:
        if(tip.fog == "none"){
            tmp <- OUwie.sim(tmp.phy, data, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, alpha=alpha, sigma.sq=sigma.sq, theta=theta, theta0=theta0, tip.fog=tip.fog, shift.point=shift.point)
        }
        if(tip.fog == "known"){
			if(!dim(data)[2]==4){
				stop("You specified tip fog should be incorporated, but this information is missing.", .call=FALSE)
			}else{
				tmp <- OUwie.sim(tmp.phy, data[,c(1,2,4)], simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, alpha=alpha, sigma.sq=sigma.sq, theta=theta, theta0=theta0, tip.fog=tip.fog, shift.point=shift.point)
			}
        }
		#OUwie.sim outputs the trait file in the order of the tree, but the trait file is likely not to be this way. So I alphabetized the input trait file above, and I do the same to the simulated trait file:
		data <- data[order(data[,1]),]
		tmp <- tmp[order(tmp[,1]),]
		#Replaces the data column with the simulated data:
        if(simmap.tree == TRUE){
            data[,3] <- tmp[,2]
        }else{
            data[,3] <- tmp[,3]
        }
		#Now run OUwie, using the measurement error if it is contained within the data, to estimate the parameters from the simulated data:
		tmp <- OUwie(phy, data, model=model, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, root.station=root.station, get.root.theta=get.root.theta, shift.point=shift.point, clade=clade, tip.fog=tip.fog, diagn=diagn, quiet=quiet, warn=warn, check.identify=FALSE, algorithm=algorithm)
		#Now bind all the relevant output together
		if(model == "BM1" | model == "BMS" | model == "OU1"){
			res <- rbind(res, c(tmp$solution[1,], tmp$solution[2,], tmp$theta[1,1]))
		}else{
			res <- rbind(res, c(tmp$solution[1,], tmp$solution[2,], tmp$theta[,1]))
		}
	}
    if(model=="BM1" | model=="OU1"){
        root.station=TRUE
    }
    if(model=="BMS"){
        root.station=FALSE
    }
	if(root.station==TRUE){
		theta.mat<-matrix(t(tmp$theta), 2, length(levels(tmp$tot.states)))
		rownames(theta.mat)<-c("estimate", "se")
		if(tmp$simmap.tree==FALSE){
			colnames(theta.mat)<- levels(tmp$tot.states)
		}
		if(tmp$simmap.tree==TRUE){
			colnames(theta.mat) <- c(colnames(tmp$phy$mapped.edge))
		}
		if(model=="BM1" | model=="OU1"){
			colnames(res) <- c(paste("alpha", levels(tmp$tot.states),sep="_"), paste("sigma.sq", levels(tmp$tot.states), sep="_"), paste("theta", "root", sep="_"))
		}else{
			colnames(res) <- c(paste("alpha", levels(tmp$tot.states),sep="_"), paste("sigma.sq", levels(tmp$tot.states), sep="_"), paste("theta", colnames(theta.mat), sep="_"))
		}
	}else{
        if(get.root.theta == TRUE){
            theta.mat<-matrix(t(tmp$theta), 2, length(levels(tmp$tot.states))+1)
            rownames(theta.mat)<-c("estimate", "se")
            if(tmp$simmap.tree==FALSE){
                colnames(theta.mat)<-c("root", levels(tmp$tot.states))
            }
            if(tmp$simmap.tree==TRUE){
                colnames(theta.mat)<-c("root", colnames(tmp$phy$mapped.edge))
            }
        }else{
            theta.mat<-matrix(t(tmp$theta), 2, length(levels(tmp$tot.states)))
            rownames(theta.mat)<-c("estimate", "se")
            if(tmp$simmap.tree==FALSE){
                colnames(theta.mat)<-c(levels(tmp$tot.states))
            }
            if(tmp$simmap.tree==TRUE){
                colnames(theta.mat)<-c(colnames(tmp$phy$mapped.edge))
            }
        }
		if(model=="BMS"){
			colnames(res) <- c(paste("alpha", levels(tmp$tot.states), sep="_"), paste("sigma.sq", levels(tmp$tot.states), sep="_"),paste("theta", "root", sep="_"))
		}else{
			colnames(res) <- c(paste("alpha", levels(tmp$tot.states), sep="_"), paste("sigma.sq", levels(tmp$tot.states), sep="_"),paste("theta", colnames(theta.mat), sep="_"))
		}
	}
	class(res) <- "OUwie.boot"
	return(res)
}



