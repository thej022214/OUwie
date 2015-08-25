\name{OUwie.fixed}
\alias{OUwie.fixed}
\title{Generalized Hansen model likelihood calculator}
\description{Allows the user to calculate the likelihood given a specified set of parameter values}
\usage{
OUwie.fixed(phy, data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"), 
simmap.tree=FALSE, scaleHeight=FALSE, root.station=TRUE, alpha=NULL, 
sigma.sq=NULL, theta=NULL, clade=NULL, mserr="none", quiet=FALSE)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes}
\item{data}{a dataframe containing species information (see Details)}
\item{model}{models to fit to comparative data (see Details).}
\item{simmap.tree}{a logical indicating whether the input tree is in SIMMAP format. The default is \code{FALSE}.} 
\item{scaleHeight}{a logical indicating whether the total tree height should be scaled to 1 (see Details). The default is \code{FALSE}.} 
\item{root.station}{a logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated.}
\item{alpha}{a numeric vector giving the values of \eqn{\alpha}{alpha} for each selective regime}
\item{sigma.sq}{a numeric vector giving the values of \eqn{\sigma^2}{sigma^2} for each selective regime} 
\item{theta}{a numeric vector giving the values of \eqn{\theta}{theta} for each selective regime}
\item{clade}{a list containing a pair of taxa whose MRCA is the clade of interest.}
\item{mserr}{designates whether a fourth column in the data matrix contains measurement error for each species value ("known"). The measurement error is assumed to be the standard error of the species mean. The default is "none".}
\item{quiet}{a logical indicating whether or not to print progress to the screen. The default is "FALSE".}
}
\details{
The input is a tree and a data file. The tree must be of class \dQuote{phylo} and must contain the ancestral selective regimes as internal node labels. The data file is a dataframe that must have column entries in the following order: [,1] species names and [,2] their current selective regime. The user specifies the  parameter values (i.e. \eqn{\alpha}{alpha}, \eqn{\sigma^2}{sigma^2}, and \eqn{\theta}{theta}).}

\value{
\code{OUwie.fixed} returns an object of class \code{OUwie.fixed}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$model}{The model being fit}
\item{$solution}{a matrix containing the maximum likelihood estimates of \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2}.}
\item{$theta}{a matrix containing the maximum likelihood estimates of \eqn{\theta}{theta}.}
\item{$tot.state}{A vector of names for the different regimes}
\item{$index.mat}{The indices of the parameters being estimated are returned. The numbers correspond to the row in the \code{eigvect} and can useful for identifying the parameters that are causing the objective function to be at a saddlepoint (see Details)}
\item{$simmap.tree}{A logical indicating whether the input phylogeny is a SIMMAP formatted tree.}
\item{$data}{User-supplied dataset}
\item{$phy}{User-supplied tree}
\item{$root.station}{A logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, was estimated}
\item{$res}{A vector of residuals from the model fit.  The residuals are ordered in the same way as the tips in the tree.}
}

\examples{
data(tworegime)

#Calculate the likelihood based on known values of
#alpha, sigma^2, and theta:
alpha=c(0.5632459,0.1726052)
sigma.sq=c(0.1064417,0.3461386)
theta=c(1.678196,0.4185894)

OUwie.fixed(tree,trait,model=c("OUMVA"), simmap.tree=FALSE, scaleHeight=FALSE,
clade=NULL, alpha=alpha,sigma.sq=sigma.sq,theta=theta)
}
\author{Jeremy M. Beaulieu}
