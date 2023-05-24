#OUwie likelihood calculator

#written by Jeremy M. Beaulieu and Brian O'Meara

#Allows the user to calculate the likelihood given a specified set of parameter values while estimating the states at internal nodes. Assumes you have estimated the paramters already using OUwie

# Idea is to take a tree, add tips at each node, add tips to the data (including figuring out regimes), and optimize the data for these made up tips. Take in OUwie object

attach.stub.taxon <- function(node, phy, tip.name=NULL) {
  if(is.null(tip.name)) {
    tip.name <- paste0("node_", node)
  }
  return(ape::bind.tree(phy, structure(list(edge = structure(c(2L, 1L), .Dim = 1:2), tip.label = tip.name, Nnode = 1L, edge.length = 0), .Names = c("edge", "tip.label", "Nnode", "edge.length"), class = "phylo"), where=node))
}

attach.stub.taxa <- function(phy) {
    root.node <- ape::Ntip(phy)+1
    start.node <- root.node+1
    n.node.not.root <- ape::Nnode(phy)-1
    for (i in sequence(n.node.not.root)) {
        phy <- attach.stub.taxon(start.node + 2*(i-1), phy, tip.name = paste0("node_", root.node+i)) #all the node numbers go up by one as we add tips
    }
    return(phy)
}

add.stub.taxon.to.data <- function(taxon, phy, data) {
    data[,1] <- as.character(data[,1])
    new.row <- data[1,]
    new.row[1,] <- NA
    new.row[1,1] <- taxon
    if(ncol(new.row)==4) { #we have measurement error
        new.row[1,4] <- 0 #but not in our internal node
    }
    taxon.id <- which(phy$tip.label==taxon)
    focal.node.id <- phy$edge[which(phy$edge[,2]==taxon.id),1]
    new.row[1,2] <- phy$node.label[focal.node.id - ape::Ntip(phy)]
    return(rbind(data, new.row))
}

add.stub.taxa.to.data <- function(phy, data) {
    taxa.to.add <- phy$tip.label[grepl("node_", phy$tip.label)]
    for (i in sequence(length(taxa.to.add))) {
        data <- add.stub.taxon.to.data(taxa.to.add[i], phy, data)
    }
    data <- data[match(phy$tip.label, data[,1]), ]
    return(data)
}

#So, OUwie data going in has a col for taxon names. Coming out of OUwie (and saved in OUwie object), it has rownames instead. So use these functions. Convert the OUwie object to have a regular OUwie formatted input
ouwie.col.to.row.names <- function(x) {
    rownames(x) <- x[,1]
    return(x[,-1])
}

ouwie.row.names.to.col <- function(x) {
    return(cbind(Genus_species=rownames(x), x))
}



anc.likelihood <- function(x, fitted.OUwie.object) {
    traits <- fitted.OUwie.object$data
    traits[grepl("node_", traits[,1]),3] <- x
    #model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"),simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE,root.station=TRUE, alpha=NULL, sigma.sq=NULL, theta=NULL, clade=NULL, mserr="none", quiet=FALSE)
    return(-1*OUwie.fixed(phy=fitted.OUwie.object$phy,data=traits, model=fitted.OUwie.object$model, simmap.tree=fitted.OUwie.object$simmap.tree, root.age=fitted.OUwie.object$root.age, scaleHeight=FALSE, root.station=fitted.OUwie.object$root.station, shift.point=fitted.OUwie.object$shift.point,  alpha=fitted.OUwie.object$solution['alpha',], sigma.sq=fitted.OUwie.object$solution['sigma.sq',], theta=fitted.OUwie.object$theta[,1], clade=NULL, mserr=ifelse(is.null(fitted.OUwie.object$mserr.est), "none", fitted.OUwie.object$mserr.est), quiet=TRUE, check.identify=FALSE, algorithm=fitted.OUwie.object$algorithm)$loglik)
}


OUwie.anc <- function(fitted.OUwie.object, opts = list("algorithm"="NLOPT_LN_BOBYQA", "maxeval"="1000", "ftol_abs"=0.001), knowledge=FALSE, multiple_starts=1){
    if(!knowledge) {
      stop("You are trying to run the function without having read all the documentation. Please do ?OUwie.anc before using this function and read ALL the documentation.")
    }
	
	fitted.OUwie.object$phy <- attach.stub.taxa(fitted.OUwie.object$phy)
    fitted.OUwie.object$data <- add.stub.taxa.to.data(fitted.OUwie.object$phy, ouwie.row.names.to.col(fitted.OUwie.object$data))
    traits <- fitted.OUwie.object$data
	n_values <- sum(grepl("node_", traits[,1]))
	recon_values <- rep(NA, n_values)
	names(recon_values) <- traits[grepl("node_", traits[,1]),1]
	print("Doing reconstructions")
	for (i in sequence(n_values)) {
		cat("Working on node ", i, " of ", n_values, "\r")
		to_delete <- names(recon_values)[-i]
		local.fitted <- fitted.OUwie.object
		local.fitted$phy <- ape::drop.tip(local.fitted$phy, to_delete)
		local.fitted$data <- local.fitted$data[-match(to_delete, local.fitted$data[,1]),]
		starting_val <- median(traits[,3], na.rm=TRUE)
		result <- nloptr(x0=starting_val, eval_f=anc.likelihood, opts=opts, fitted.OUwie.object=local.fitted)
		recon_values[i] <- result$solution
		for (j in sequence(multiple_starts-1)) {
			algorithm_original <- opts$algorithm
			opts$algorithm <- sample(c("NLOPT_LN_SBPLX", "NLOPT_LN_BOBYQA"), 1)
			new_starting_val <- rnorm(1, mean=sample(c(starting_val, recon_values[i]), 1), sd=0.1*abs(recon_values[i]))
			new_result <- nloptr(x0=new_starting_val, eval_f=anc.likelihood, opts=opts, fitted.OUwie.object=local.fitted)
			opts$algorithm <- algorithm_original
			if(new_result$objective < result$objective) {
				result <- new_result
				recon_values[i] <- result$solution
			}
		}
	}
	print("Done reconstructions")
	traits[grepl("node_", traits[,1]),3] <- recon_values
 	fitted.OUwie.object$data <- ouwie.col.to.row.names(traits)
    quantitative.trait <- fitted.OUwie.object$data[,2]
    names(quantitative.trait) <- rownames(fitted.OUwie.object$data)
    user.recons <- rep(NA, ape::Nnode(fitted.OUwie.object$phy))
    for(i in sequence(ape::Nnode(fitted.OUwie.object$phy))) {
        if(i==1) { #we're at the root
            user.recons[1] <- fitted.OUwie.object$theta[min(nrow(fitted.OUwie.object$theta), fitted.OUwie.object$phy$node.label[1]),1] #so get the estimate for the regime of there are more than one, otherwise, BM
        } else {
            node.index <- i+ape::Ntip(fitted.OUwie.object$phy)-sum(grepl("node_", fitted.OUwie.object$phy$tip.label))
            # match up here
            user.recons[i] <- fitted.OUwie.object$data[paste0("node_", node.index),2]
        }
    }
    fitted.OUwie.object$NodeRecon <- user.recons
    class(fitted.OUwie.object) <- c("OUwie.anc", "OUwie")
    return(fitted.OUwie.object)
}



OUwie.anc.older.simultaneous <- function(fitted.OUwie.object, opts = list("algorithm"="NLOPT_LN_BOBYQA", "maxeval"="1000", "ftol_abs"=0.001), knowledge=FALSE, multiple_starts=1){
    if(!knowledge) {
      stop("You are trying to run the function without having read all the documentation. Please do ?OUwie.anc before using this function and read ALL the documentation.")
    }
	
	fitted.OUwie.object$phy <- attach.stub.taxa(fitted.OUwie.object$phy)
    fitted.OUwie.object$data <- add.stub.taxa.to.data(fitted.OUwie.object$phy, ouwie.row.names.to.col(fitted.OUwie.object$data))
    traits <- fitted.OUwie.object$data
	n_values <- sum(grepl("node_", traits[,1]))
	print(paste0("Starting optimization 1 of ", multiple_starts))
	median_values <- rep(median(traits[,3], na.rm=TRUE), n_values)
	print(median_values)
    result <- nloptr(x0=median_values, eval_f=anc.likelihood, opts=opts, fitted.OUwie.object=fitted.OUwie.object)
	if(multiple_starts>1) {
		print(paste0("Found reconstruction with likelihood: ", result$objective))
		print(result$solution)
		for(i in sequence(multiple_starts-1)) {
			opts$algorithm <- sample(c("NLOPT_LN_SBPLX", "NLOPT_LN_BOBYQA"), 1)
			print(paste0("Starting optimization ", i+1, " of ", multiple_starts))
			x0_start <- rnorm(n=n_values, mean=ifelse(sample(c(TRUE, FALSE), n_values, replace=TRUE, prob=c(0.8, 0.2)), median_values, result$solution), sd=runif(1, min=0, max=sqrt(diff(range(traits[,3], na.rm=TRUE)))))
			print(x0_start)
			result_restart <- nloptr(x0=x0_start, eval_f=anc.likelihood, opts=opts, fitted.OUwie.object=fitted.OUwie.object)
			print(paste0("Found reconstruction with likelihood: ", result_restart$objective))
			print(result_restart$solution)
			if(result_restart$objective < result$objective) {
				print(paste0("Found better reconstruction: ", result_restart$objective, " vs ", result$objective))
				result <- result_restart
			}
		}	
	}
    traits[grepl("node_", traits[,1]),3] <- result$solution
    fitted.OUwie.object$data <- ouwie.col.to.row.names(traits)
    quantitative.trait <- fitted.OUwie.object$data[,2]
    names(quantitative.trait) <- rownames(fitted.OUwie.object$data)
    user.recons <- rep(NA, ape::Nnode(fitted.OUwie.object$phy))
    for(i in sequence(ape::Nnode(fitted.OUwie.object$phy))) {
        if(i==1) { #we're at the root
            user.recons[1] <- fitted.OUwie.object$theta[min(nrow(fitted.OUwie.object$theta), fitted.OUwie.object$phy$node.label[1]),1] #so get the estimate for the regime of there are more than one, otherwise, BM
        } else {
            node.index <- i+ape::Ntip(fitted.OUwie.object$phy)-sum(grepl("node_", fitted.OUwie.object$phy$tip.label))
            # match up here
            user.recons[i] <- fitted.OUwie.object$data[paste0("node_", node.index),2]
        }
    }
    fitted.OUwie.object$NodeRecon <- user.recons
    class(fitted.OUwie.object) <- c("OUwie.anc", "OUwie")
    return(fitted.OUwie.object)
}

plot.OUwie.anc <- function(x, ...) {
    quantitative.trait <- x$data[,2]
    names(quantitative.trait) <- rownames(x$data)
    # user.recons <- rep(NA, ape::Nnode(x$phy))
    # for(i in sequence(ape::Nnode(x$phy))) {
    #     if(i==1) { #we're at the root
    #         user.recons[1] <- x$theta[min(nrow(x$theta), x$phy$node.label[1]),1] #so get the estimate for the regime of there are more than one, otherwise, BM
    #     } else {
    #         node.index <- i+ape::Ntip(x$phy)-sum(grepl("node_", x$phy$tip.label))
    #         # match up here
    #         user.recons[i] <- x$data[paste0("node_", node.index),2]
    #     }
    # }
    x$phy <- ape::drop.tip(x$phy, x$phy$tip.label[grepl("node_",x$phy$tip.label)])
    pruned <- geiger::treedata(x$phy, quantitative.trait, warnings=FALSE, sort=TRUE)
    phytools::contMap(pruned$phy, pruned$data[,1], method="user", anc.states=x$NodeRecon, ...)
}
