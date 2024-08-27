\name{OUwie.dredge}
\alias{OUwie.dredge}
\title{Generalized Detection of shifts in OU process}
\description{Allows the hypothesis free detection of shifts in the OU process. The number and location of shifts is estimated using a user-specified information criterion.}
\usage{
OUwie.dredge(phy, data, criterion=c("AIC", "AICc", "BIC", "mBIC"), shift.max=3, 
sigma.sq.max.k=3, alpha.max.k=3, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, 
shift.point=0.5, tip.fog="none", algorithm=c("invert", "three.point"), lb=NULL, ub=NULL, 
opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000", 
"ftol_rel"=.Machine$double.eps^0.5), verbose=FALSE)
}  
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes}
\item{data}{a dataframe containing two columns, taxon names in the first column, and species trait information in the second column.}
\item{criterion}{information criterion to use for shift detection.}
\item{shift.max}{maximum allowed number of shifts.} 
\item{sigma.sq.max.k}{maximum allowed number of sigma.sq parameters.} 
\item{alpha.max.k}{maximum allowed number of alpha parameters.} 
\item{root.age}{indicates the age of the tree. This is to be used in cases where the "tips" are not contemporary, such as in cases for fossil trees. Default is \code{NULL} meaning latest tip is modern day.}
\item{scaleHeight}{a logical indicating whether the total tree height should be scaled to 1. The default is \code{FALSE}.} 
\item{root.station}{a logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated.}
\item{shift.point}{the point along a branch where a regime change is assumed to have occurred (if \code{SIMMAP=FALSE}. The default is set to 0.5, or halfway along a branch.}
\item{tip.fog}{designates whether a fourth column in the data matrix contains tip fog for each species value ("known"). The tip fog is assumed to be the standard error of the species mean. The default is "none".}
\item{algorithm}{designates whether the standard matrix inversion ('invert') or the faster 'three-point' algorithm of Ho and Ane (2013) should be used.}
\item{lb}{if algorithm == "invert" a single value indicating the lower bound for the parameter values. if algorithm == "three.point", a vector of length 3 with position 1 as alpha lower bound, position 2 as sigma.sq's lower bound, position 3 as theta's lower bound. when set to NULL it will be a default value of 1e-9. Note that even if the model you're using doesn't include alpha (e.g. BM1), it must be included in this vector, but it will not be used to set any bounds.}
\item{ub}{if algorithm == "invert" a single value indicating the upper bound for the parameter values. if algorithm == "three.point", a vector of length 3 with position 1 as alpha upper bound, position 2 as sigma.sq's upper bound, position 3 as theta's upper bound. when set to NULL it will be a default value of 100. Note that even if the model you're using doesn't include alpha (e.g. BM1), it must be included in this vector, but it will not be used to set any bounds.}
\item{opts}{a list of options to pass to nloptr for the optimization: useful to adjust for faster, coarser searches}
\item{verbose}{a logical indicating whether to print incremental steps}
}
\details{
This is an expanded version of the shift point model of Ho and Ane (2014). This is currently being tested, but as of now we strongly recommend using the mBIC criterion when searching for shifts.
}
\value{
\code{OUwie.dredge} returns an object of class \code{OUwie.dredge}. This is a list with elements:
\item{$loglik}{the maximum log-likelihood.}
\item{$criterion}{the information criterion to use for shift detection.}
\item{$criterion.score}{the information criterion score used for shift detection.}
\item{$shift.model}{The shift model estimated from the data.}
\item{$solution}{a matrix containing the maximum likelihood estimates of \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2}.}
\item{$tip.fog.est}{indicates value of the measurement error if it was estimated from the data.}
\item{$theta}{a matrix containing the maximum likelihood estimates of \eqn{\theta}{theta}.}
\item{$tot.states}{A vector of names for the different regimes}
\item{$index.mat}{The indices of the parameters being estimated are returned. The numbers correspond to the row in the \code{eigvect} and can useful for identifying the parameters that are causing the objective function to be at a saddlepoint (see Details)}
\item{$simmap.tree}{A logical indicating whether the input phylogeny is a SIMMAP formatted tree.}
\item{$root.age}{The user-supplied age at the root of the tree.}
\item{$scaleHeight}{Indicates whether the tree was constrained to a total height of 1.}
\item{$shift.point}{The user-specified portion of the branch where a regime shift occurs.}
\item{$opts}{Settings used for optimization routine.}
\item{$data}{The shift model dataset, which includes regime painting for each tip.}
\item{$phy}{The shift model painted phylogeny.}
\item{$root.station}{A logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, was estimated}
\item{$starting.vals}{the starting values used for the parameter search.}
\item{$regime.weights}{A table containing parameter estimates and the weights for time spent in each regime for each tip.}
}
\references{
Ho, L.S.T., and C. Ane. 2014. Intrinsic inference difficulties for trait evolution with Ornstein-Uhlenbeck models. Methods in Ecology and Evolution, 5: 1133-1146.
}
\author{Jeremy M. Beaulieu}
