\name{OUwie.slice}
\alias{OUwie.slice}
\title{Generalized Hansen models with time slices}
\description{Fits generalized Ornstein-Uhlenbeck-based Hansen models of continuous characters before and after time slices.}
\usage{
OUwie.slice(phy, data, model=c("BMS","OUM","OUMV","OUMA","OUMVA"),
 timeslices=c(NA), scaleHeight=FALSE, root.station=TRUE, mserr="none",
 diagn=FALSE, quiet=FALSE, warn=TRUE)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data matrix containing species information (see Details).}
\item{model}{models to fit to comparative data (see Details).}
\item{timeslices}{specifies the value and number fixed timeslices, timeslices to be estimated, or both (see Details).}
\item{scaleHeight}{a logical indicating whether the total tree height should be scaled to 1 (see Details). The default is \code{FALSE}.} 
\item{root.station}{a logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated (see Details).} 
\item{mserr}{designates whether a fourth column in the data matrix contains measurement error for each species value ("known"). The measurement error is assumed to be the standard error of the species mean. The default is "none".}
\item{diagn}{a logical indicating whether the full diagnostic analysis should be carried out. The default is \code{FALSE}.}
\item{quiet}{a logical indicating whether progress should be written to the screen. The default is \code{FALSE}.}
\item{warn}{a logical indicating whether a warning should be printed if the number of parameters exceeds ntips/10. The default is \code{TRUE}.}
}
\details{
This function fits various likelihood models for continuous characters evolving under discrete selective regimes that defined by a time slice (i.e., before and after K-P event). Timeslices can be fixed, estimated from the data, or some combination of both. Here it is assumed that the present is T=0 and the root is max(branching.times(phy)). Also note that this function is still under development and so far it seems that when specifying an OU model there is a narrow range of conditions in which meaningful parameter estimates are obtained. Thus, use this function at your own risk.

Also note, when specifying the BMS model be mindful of the root.station flag. When root.station=FALSE, the non-censored model of O'Meara et al. 2006 is invoked (i.e., a single regime at the root is estimated), and when root.station==TRUE the group mean model of Thomas et al. 2006 (i.e., the number of means equals the number of regimes). The latter case appears to be a strange special case of OU, in that it behaves similarly to the OUMV model, but without selection. I would say that this is more consistent with the censored test of O'Meara et al. (2006), as opposed to having any real connection to OU. In any case, more work is clearly needed to understand the behavior of the group means model, and therefore, I recommend setting root.station=FALSE in the BMS case.
}
\value{
\code{OUwie.slice} returns an object of class \code{OUwie.slice}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$model}{The model being fit}
\item{$solution}{a matrix containing the maximum likelihood estimates of \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2}.}
\item{$theta}{a matrix containing the maximum likelihood estimates of \eqn{\theta}{theta} and it standard error.}
\item{$solution.se}{a matrix containing the approximate standard errors of \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2}. The standard error is calculated as the diagonal of the inverse of the Hessian matrix.}
\item{$timeslices}{a vector of timeslices either based on fixed age specified by the user, estimated from the data, or both.}
\item{$tot.state}{A vector of names for the different regimes}
\item{$index.mat}{The indices of the parameters being estimated are returned. The numbers correspond to the row in the \code{eigvect} and can useful for identifying the parameters that are causing the objective function to be at a saddlepoint (see Details)}
\item{$simmap.tree}{A logical indicating whether the input phylogeny is a SIMMAP formatted tree.}
\item{$opts}{Internal settings of the likelihood search}
\item{$data}{User-supplied dataset}
\item{$phy}{User-supplied tree}
\item{$root.station}{A logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, was estimated}
\item{$lb}{The lower bound set}
\item{$ub}{The upper bound set}
\item{$iterations}{Number of iterations of the likelihood search that were executed}
\item{$mserr.est}{The estimated measurement error if mserr="est". Otherwise, the value is NULL.}
\item{$res}{A vector of residuals from the model fit.  The residuals are ordered in the same way as the tips in the tree.}
\item{$eigval}{The eigenvalues from the decomposition of the Hessian of the likelihood function. If any \code{eigval<0} then one or more parameters were not optimized during the likelihood search (see Details)}
\item{$eigvect}{The eigenvectors from the decomposition of the Hessian of the likelihood function is returned (see Details)}
}
\examples{
data(tworegime)

##Here we want a fixed slice at T=2, assuming the present is T=0:
#library(phytools)
#max.height <- max(nodeHeights(tree))
#timeslices <- max.height - 2
#timeslices <- c(0,timeslices)
#phy.sliced<-make.era.map(tree,timeslices)
#leg<-c("blue3","red3")
#names(leg)<-c(1,2)
#plotSimmap(phy.sliced,leg, pts=FALSE, ftype="off", lwd=1)

##Now fit an BMS model with a single fixed timeslice at time=2:
#ppBM<-OUwie.slice(tree,trait[,c(1,3)],model=c("BMS"), root.station=TRUE, timeslices=c(2))

##Fit an OU model with a single fixed timeslice:
#ppOUM<-OUwie.slice(tree,trait[,c(1,3)],model=c("OUM"), root.station=TRUE, timeslices=c(2))

##Fit an BMS model with an unknown timeslice:
#ppBM<-OUwie.slice(tree,trait[,c(1,3)],model=c("BMS"), root.station=TRUE, timeslices=c(NA))

##Fit an BMS model with an unknown and a fixed timeslice:
#ppBM<-OUwie.slice(tree,trait[,c(1,3)],model=c("BMS"), root.station=TRUE, timeslices=c(NA,2))

}
\references{
Beaulieu J.M., and O'Meara B.C. In prep.

Beaulieu J.M., Jhwueng D.C., Boettiger C., and O'Meara B.C. 2012. Modeling stabilizing selection: Expanding the Ornstein-Uhlenbeck model of adaptive evolution. Evolution 66:2369-2383.

O'Meara B.C., Ane C., Sanderson P.C., Wainwright P.C. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60:922-933.

Butler M.A., King A.A. 2004. Phylogenetic comparative analysis: A modeling approach for adaptive evolution. American Naturalist 164:683-695.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
