\name{OUwie.contour}
\alias{OUwie.contour}
\title{Contour plot of likelihood surface}
\description{Plots likelihood surface for pairs of free parameters for generalized Ornstein-Uhlenbeck-based Hansen models of continuous characters evolving under discrete selective regimes.}
\usage{
OUwie.contour(phy, data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"), 
simmap.tree=FALSE, scaleHeight=FALSE,root.station=TRUE, lb=0.000001, ub=1000, 
focal.param=NULL, clade=NULL, mserr="none", nrep=1000, sd.mult=3, 
levels=c(0.5,1,1.5,2), likelihood.boundary=Inf, lwd=2, \dots)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data matrix containing species information (see Details).}
\item{model}{models to fit to comparative data (see Details).}
\item{simmap.tree}{a logical indicating whether the input tree is in SIMMAP format. The default is \code{FALSE}.} 
\item{scaleHeight}{a logical indicating whether the total tree height should be scaled to 1 (see Details). The default is \code{FALSE}.} 
\item{root.station}{a logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated (see Details).}
\item{lb}{lower bound for the likelihood search. The default is \code{lb=0.000001}. Note that this value must be greater than zero.}
\item{ub}{upper bound for the likelihood search. The default is \code{ub=1000}.}
\item{focal.param}{a character vector with the pair of parameters (sigma.sq and/or alpha) to plot (see Details).}
\item{clade}{a list containing a pair of taxa whose MRCA is the clade of interest (see Details).}
\item{mserr}{designates whether a fourth column in the data matrix contains measurement error for each species value ("known"). The measurement error is assumed to be the standard error of the species mean. The default is "none".}
\item{nrep}{the number of points to use for the likelihood surface.}
\item{sd.mult}{the sd for the normal distribution to pull random points from is sd.mult * se for each parameter.}
\item{levels}{the levels at which to draw contour lines, measured as lnL units away from the best values.}
\item{likelihood.boundary}{default = Inf}
\item{lwd}{a graphical control for the width of the line}
\item{\dots}{Additional arguments to be passed to the plot device} 
}
\details{
This function creates a plot of the likelihood surface for a pair of parameters, letting the other parameters find their own optima. It can be very slow, as it involves optimization \code{nrep} times (though with two fewer parameters than with the chosen model, as the focal parameter values are fixed). It uses half its points to sample near the MLE and half to uniformly sample across a range. However, it is important to do as on some trees and with some models, certain parameters may correlate with each other.

The pair of parameters to examine is passed by focal.param. For example, to do a plot of sigma.sq from the first regime and alpha from the second regime, one would pass 

\code{focal.param = c( "sigma.sq_1", "alpha_2")}.

This returns a data.frame with the first two columns being the values of the points examined and the third column the loglik of those points. The last row contains the MLE.

IMPORTANT: this currently written to only do sigma.sq and/or alpha, not theta.
IMPORTANT II: this function requires the package akima.
}
\value{
\item{finalResults}{the parameter values and loglik}
}
\references{
Beaulieu J.M., Jhwueng D.C., Boettiger, C., and O'Meara B.C. 2012. Modeling stabilizing selection: Expanding the Ornstein-Uhlenbeck model of adaptive evolution. Evolution 66:2369-2383.

O'Meara B.C., Ane C., Sanderson P.C., Wainwright P.C. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60:922-933.

Butler M.A., King A.A. 2004. Phylogenetic comparative analysis: A modeling approach for adaptive evolution. American Naturalist 164:683-695.
}
\author{Brian C. O'Meara and Jeremy M. Beaulieu}
\keyword{models}

