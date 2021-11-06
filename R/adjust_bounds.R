fix_upper <- function(ub, x0) {
	if(any(ub<x0)) {
		warning("Some upper bounds are less than initial values, adjusting")
		x0[which(ub<x0)] <- ub[which(ub<x0)]
	}	
	return(ub)
}

fix_lower <- function(lb, x0) {
	if(any(lb>x0)) {
		warning("Some lower bounds are greater than initial values, adjusting")
		x0[which(lb>x0)] <- lb[which(lb>x0)]
	}	
	return(lb)
}