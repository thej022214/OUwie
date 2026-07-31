#Weight matrix generator taken from Butler and King (2004) and modified to allow multiple alpha parameters
 
#written by Jeremy M. Beaulieu

weight.mat <- function(phy, edges, Rate.mat, root.state, simmap.tree = FALSE, root.age = NULL, scaleHeight = FALSE, assume.station = TRUE, shift.point = 0.5){
	age.table <- MakeAgeTable(phy, root.age = root.age)
	Tmax <- max(age.table)
	n <- max(phy$edge[, 1])
	ntips <- length(phy$tip.label)
	
	if(is.null(root.state)) {
		root.state <- which(edges[dim(edges)[1], ] == 1) - 5
		edges <- edges[-1 * dim(edges)[1], ]
	}

	if(simmap.tree == TRUE) {
		k <- length(colnames(phy$mapped.edge))
	}else{
		mm <- dim(edges)
		k <- length(6:mm[2])
	}
	alpha <- Rate.mat[1, ]

	root_node <- setdiff(unique(edges[, 2]), unique(edges[, 3]))[1]
	max_node_id <- max(c(edges[, 2], edges[, 3]))
	Ato <- numeric(max_node_id)
	Ato[root_node] <- 0

	nodevar.root.tot <- rep(0, max(edges[, 3]))
	nodevar.k <- rep(0, max(edges[, 3]))

	W <- matrix(0, ntips, k)

	#Ato[ancestor] must be filled before the edge below it is visited. edges and
	#phy$maps share an index, so we visit that index in cladewise order. a NULL
	#root.state drops a row from edges above, so bound the traversal by what is left.
	edge.order <- ape::reorder.phylo(phy, "cladewise", index.only = TRUE)
	edge.order <- edge.order[edge.order <= nrow(edges)]

	#the regime and duration of each epoch, the Ato accumulation and nodevar.root.tot
	#are all the same for every regime j, so they are built once here instead of being
	#recomputed k times inside the loop below
	seg.regime <- vector("list", nrow(edges))
	seg.dt <- vector("list", nrow(edges))
	Astart <- numeric(nrow(edges))
	n.cov.root.tot <- matrix(0, n, 1)

	for(i in edge.order){
		anc <- edges[i, 2]
		desc <- edges[i, 3]

		#cumulative A at start of the edge
		Acur <- Ato[anc]
		Astart[i] <- Acur
		nodevar.root.tot[i] <- 0

		if(simmap.tree == TRUE){
			if(scaleHeight == TRUE) {
				currentmap <- phy$maps[[i]] / Tmax
			}else{
				currentmap <- phy$maps[[i]]
			}
			regimes <- numeric(length(currentmap))
			dts <- numeric(length(currentmap))
			for(regimeindex in 1:length(currentmap)){
				dts[regimeindex] <- as.numeric(currentmap[regimeindex])
				regimes[regimeindex] <- which(colnames(phy$mapped.edge) == names(currentmap)[regimeindex])
			}
		}else{
			oldtime <- edges[i, 4]
			newtime <- edges[i, 5]
			if(anc %in% edges[, 3]){
				start <- which(edges[, 3] == anc)
				oldregime <- which(edges[start, 6:(k + 5)] == 1)
			}else{
				oldregime <- root.state
			}
			newregime <- which(edges[i, 6:(k + 5)] == 1)
			if(oldregime == newregime){
				regimes <- newregime
				dts <- newtime - oldtime
			}else{
				shifttime <- newtime - ((newtime - oldtime) * shift.point)
				# epoch 1 - oldregime, epoch 2 - newregime
				regimes <- c(oldregime, newregime)
				dts <- c(shifttime - oldtime, newtime - shifttime)
			}
		}

		for(regimeindex in seq_along(dts)){
			a <- alpha[regimes[regimeindex]]
			nodevar.root.tot[i] <- nodevar.root.tot[i] - a * dts[regimeindex]
			Acur <- Acur + a * dts[regimeindex]
		}

		seg.regime[[i]] <- regimes
		seg.dt[[i]] <- dts
		n.cov.root.tot[desc, ] <- nodevar.root.tot[i]

		Ato[desc] <- Acur
	}

	w.root.tot <- mat.gen.diag(phy, n.cov.root.tot)
	w_root <- exp(w.root.tot)

	#only nodevar.k depends on the regime being weighted
	for(j in 1:k){
		n.cov.k <- matrix(0, n, 1)
		for(i in edge.order){
			Acur <- Astart[i]
			nodevar.k[i] <- 0
			regimes <- seg.regime[[i]]
			dts <- seg.dt[[i]]
			for(regimeindex in seq_along(dts)){
				a <- alpha[regimes[regimeindex]]
				if(regimes[regimeindex] == j){
					nodevar.k[i] <- nodevar.k[i] + (exp(Acur + a * dts[regimeindex]) - exp(Acur))
				}
				Acur <- Acur + a * dts[regimeindex]
			}
			n.cov.k[edges[i, 3], ] <- nodevar.k[i]
		}
		W[, j] <- w_root * mat.gen.diag(phy, n.cov.k)
	}

	if (assume.station == TRUE) {
		W[, root.state] <- W[, root.state] + w_root
	}else{
		W <- cbind(w_root, W)
	}
	#Restandardizes W so that the rows sum to 1 -- Generalized. Will reduce to the simpler model if assuming 1 alpha parameter, but when alpha varies by regime they will sum to 1 (though proportionally should be ok).
	W <- W / rowSums(W)
	W
}


weight.mat.old <-function(phy, edges, Rate.mat, root.state, simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE, assume.station=TRUE, shift.point=0.5){
	
	age.table <- MakeAgeTable(phy, root.age=root.age)
	Tmax <- max(age.table)
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
					currentmap <- phy$maps[[i]]/Tmax
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

