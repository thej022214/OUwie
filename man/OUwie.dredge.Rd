\name{OUwie.dredge}
\alias{OUwie.dredge}
\title{Generalized Detection of shifts in OU process}
\description{Allows the hypothesis free detection of shifts in the OU process. The number and location of shifts is estimated using a user-specified information criterion.}
\usage{
OUwie.dredge(phy, data, criterion=c("AIC", "AICc", "BIC", "mBIC"), shift.max=3, 
sigma.sq.max.k=3, alpha.max.k=3, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, 
shift.point=0.5, mserr="none", opts = list("algorithm"="NLOPT_LN_SBPLX",
 "maxeval"="1000", "ftol_rel"=.Machine$double.eps^0.5))
}  
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes}
\item{data}{a dataframe containing species information (see Details)}
\item{criterion}{information criterion to use for shift detection (see Details).}
\item{shift.max}{maximum allowed number of shifts.} 
\item{sigma.sq.max.k}{maximum allowed number of sigma.sq parameters.} 
\item{alpha.max.k}{maximum allowed number of alpha parameters.} 
\item{root.age}{indicates the age of the tree. This is to be used in cases where the "tips" are not contemporary, such as in cases for fossil trees. Default is \code{NULL} meaning latest tip is modern day.}
\item{scaleHeight}{a logical indicating whether the total tree height should be scaled to 1 (see Details). The default is \code{FALSE}.} 
\item{root.station}{a logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated.}
\item{shift.point}{the point along a branch where a regime change is assumed to have occurred (if \code{SIMMAP=FALSE}. The default is set to 0.5, or halfway along a branch.}
\item{mserr}{designates whether a fourth column in the data matrix contains measurement error for each species value ("known"). The measurement error is assumed to be the standard error of the species mean. The default is "none".}
\item{opts}{a list of options to pass to nloptr for the optimization: useful to adjust for faster, coarser searches}
}
\value{
\code{OUwie.fixed} returns an object of class \code{OUwie.fixed}. This is a list with elements:
\item{$loglik}{the maximum log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$model}{The model being fit}
\item{$param.count}{The number of parameters counted in the model.}
\item{$solution}{a matrix containing the maximum likelihood estimates of \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2}.}
\item{$theta}{a matrix containing the maximum likelihood estimates of \eqn{\theta}{theta}.}
\item{$tot.state}{A vector of names for the different regimes}
\item{$index.mat}{The indices of the parameters being estimated are returned. The numbers correspond to the row in the \code{eigvect} and can useful for identifying the parameters that are causing the objective function to be at a saddlepoint (see Details)}
\item{$simmap.tree}{A logical indicating whether the input phylogeny is a SIMMAP formatted tree.}
\item{$root.age}{The user-supplied age at the root of the tree.}
\item{$data}{User-supplied dataset}
\item{$phy}{User-supplied tree}
\item{$root.station}{A logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, was estimated}
\item{$res}{A vector of residuals from the model fit.  The residuals are ordered in the same way as the tips in the tree.}
}
\author{Jeremy M. Beaulieu}
