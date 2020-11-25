\name{OUwie.fixed}
\alias{OUwie.fixed}
\title{Generalized Hansen model likelihood calculator}
\description{Allows the user to calculate the likelihood given a specified set of parameter values}
\usage{
OUwie.fixed(phy, data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"), 
simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, 
get.root.theta=FALSE, shift.point=0.5, alpha=NULL, sigma.sq=NULL, theta=NULL, 
clade=NULL, mserr="none", check.identify=TRUE, algorithm=c("invert", "three.point"), 
tip.paths=NULL, quiet=FALSE)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes}
\item{data}{a dataframe containing species information (see Details)}
\item{model}{models to fit to comparative data (see Details).}
\item{simmap.tree}{a logical indicating whether the input tree is in SIMMAP format. The default is \code{FALSE}.} 
\item{root.age}{indicates the age of the tree. This is to be used in cases where the "tips" are not contemporary, such as in cases for fossil trees. Default is \code{NULL} meaning latest tip is modern day.}
\item{scaleHeight}{a logical indicating whether the total tree height should be scaled to 1 (see Details). The default is \code{FALSE}.} 
\item{root.station}{a logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated.}
\item{get.root.theta}{a logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated (see Details).}
\item{shift.point}{the point along a branch where a regime change is assumed to have occurred (if \code{SIMMAP=FALSE}. The default is set to 0.5, or halfway along a branch.}
\item{alpha}{a numeric vector giving the values of \eqn{\alpha}{alpha} for each selective regime}
\item{sigma.sq}{a numeric vector giving the values of \eqn{\sigma^2}{sigma^2} for each selective regime} 
\item{theta}{a numeric vector giving the values of \eqn{\theta}{theta} for each selective regime}
\item{clade}{a list containing a pair of taxa whose MRCA is the clade of interest.}
\item{mserr}{designates whether a fourth column in the data matrix contains measurement error for each species value ("known"). The measurement error is assumed to be the standard error of the species mean. The default is "none".}
\item{check.identify}{a logical indicating whether to check that the user-supplied regime paintings will produce identifiable theta estimates. The default is \code{TRUE}.}
\item{algorithm}{designates whether matrix inversion (slow) or the three-point algorithm (fast) of Ho and Ane (2013) should be used.}
\item{tip.paths}{an optional list that can be provided by the user where each element from 1:nTips is the path from tip to root by labeled node number. The default is NULL and this option is used for internal speedups.}
\item{quiet}{a logical indicating whether or not to print progress to the screen. The default is "FALSE".}
}
\details{
The input is a tree and a data file. The tree must be of class \dQuote{phylo} and must contain the ancestral selective regimes as internal node labels. The data file is a data.frame that must have column entries in the following order: [,1] species names and [,2] their current selective regime. The user specifies the  parameter values (i.e. \eqn{\alpha}{alpha}, \eqn{\sigma^2}{sigma^2}, and \eqn{\theta}{theta}).
}
\value{
\code{OUwie.fixed} returns an object of class \code{OUwie.fixed}. This is a list with elements:
\item{$loglik}{the maximum log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$BIC}{Schwartz information criterion.}
\item{$model}{The model being fit}
\item{$param.count}{The number of parameters counted in the model.}
\item{$solution}{a matrix containing the maximum likelihood estimates of \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2}.}
\item{$theta}{a matrix containing the maximum likelihood estimates of \eqn{\theta}{theta}.}
\item{$tot.state}{A vector of names for the different regimes}
\item{$index.mat}{The indices of the parameters being estimated are returned. The numbers correspond to the row in the \code{eigvect} and can useful for identifying the parameters that are causing the objective function to be at a saddlepoint (see Details)}
\item{$simmap.tree}{A logical indicating whether the input phylogeny is a SIMMAP formatted tree.}
\item{$root.age}{The user-supplied age at the root of the tree.}
\item{$shift.point}{The user-specified portion of the branch where a regime shift occurs.}
\item{$data}{User-supplied dataset}
\item{$phy}{User-supplied tree}
\item{$root.station}{A logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, was estimated}
\item{$scaleHeight}{Indicates whether the tree was constrained to a total height of 1.}
\item{$get.root.theta}{Indicates whether the root.theta was included in the model.}
\item{$regime.weights}{A table containing parameter estimates and the weights for time spent in each regime for each tip.}
\item{$algorithm}{The algorithm used to estimate parameters.}
}

\examples{
\donttest{
data(tworegime)

#Calculate the likelihood based on known values of
#alpha, sigma^2, and theta:
alpha=c(0.5632459,0.1726052)
sigma.sq=c(0.1064417,0.3461386)
theta=c(1.678196,0.4185894)

OUwie.fixed(tree,trait,model=c("OUMVA"), simmap.tree=FALSE, scaleHeight=FALSE,
clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta, algorithm="three.point")
}
}
\author{Jeremy M. Beaulieu}
