# Format a tree and data for OUwie. By Brian O'Meara

OUwie.format <- function(phy, tip.regimes=NULL, tip.data=NULL, tip.fog=NULL, tip.fog.percentage=NULL, verbose=TRUE){
	traits <- data.frame(Taxa=phy$tip.label, Regime=1, Trait=NA, TipFog=0, stringsAsFactors=FALSE)
	if(!is.null(tip.regimes)) {
		traits$Regime <- NA # so that lack of matches don't get assigned to regime 1
		if(is.vector(tip.regimes)) {
			traits$Regime <- tip.regimes[match(traits$Taxa, names(tip.regimes))]
		} else {
			if(ncol(tip.regimes)!= 1) {
				stop("This assumes tip.regimes is either a named vector or a data.frame with one column with taxon names as rownames")	
			}
			traits$Regime <- tip.regimes[match(traits$Taxa, rownames(tip.regimes)),]
		}
		if(any(is.na(traits$Regime))) {
			warning(paste0("Some taxa (", sum(is.na(traits$Regime))," total) in the tree are not in the tip.regimes input")	)
		}
	}
	if(!is.null(tip.data)) {
		traits$Trait <- NA # so that lack of matches don't get assigned a state
		if(is.vector(tip.data)) {
			traits$Trait <- tip.data[match(traits$Taxa, names(tip.data))]
		} else {
			if(ncol(tip.data)!= 1) {
				stop("This assumes tip.data is either a named vector or a data.frame with one column with taxon names as rownames")	
			}
			traits$Trait <- tip.data[match(traits$Taxa, rownames(tip.data)),]
		}
		if(any(is.na(traits$Trait))) {
			warning(paste0("Some taxa (", sum(is.na(traits$Trait))," total) in the tree are not in the tip.data input"))
		}
	}
	if(!is.null(tip.fog)) {
		traits$TipFog <- NA # so that lack of matches don't get assigned zero
		if(is.vector(tip.fog)) {
			if(length(tip.fog)==1) {
				traits$TipFog <- tip.fog
			} else {
				traits$TipFog <- tip.fog[match(traits$Taxa, names(tip.fog))]
			}
		} else {
			if(ncol(tip.fog)!= 1) {
				stop("This assumes tip.mfog is either a named vector or a data.frame with one column with taxon names as rownames")
			}
			traits$TipFog <- tip.fog[match(traits$Taxa, rownames(tip.fog)),]
		}
		if(any(is.na(traits$TipFog))) {
			warning(paste0("Some taxa (", sum(is.na(traits$TipFog))," total) in the tree are not in the tip.fog input"))
		}
	} else if (!is.null(tip.fog.percentage)) {
		traits$TipFog <- traits$Trait * tip.fog.percentage/100
	} else {
		if(verbose) {
			warning("No tip.fog input, so assuming zero tip fog for all taxa. This is almost surely untrue in reality and is a strong bias")
		}
	}
	if(is.null(phy$node.label)) {
		phy$node.label <- rep(1, ape::Nnode(phy))	
		if(length(unique(traits$Regime[!is.na(traits$Regime)]))>1) { # we have some variation, so time to do anc recon
			if(verbose) {
				print("Doing ancestral state reconstruction of regimes, but you may want to use hOUwie instead")	
			}
			corHMM_result <- corHMM::corHMM(phy=phy, data=traits[,c("Taxa","Regime")], rate.cat=1, model='ER', node.states='joint')
			phy$node.label <- corHMM_result$states
		}
	}
	return(list(tree=phy, data=traits))
}
