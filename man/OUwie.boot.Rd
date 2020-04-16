\name{OUwie.boot}
\alias{OUwie.boot}
\title{Parametric bootstrap function}
\description{A function that performs a parametric bootstrap for a set of user-specified model parameters}
\usage{
OUwie.boot(phy, data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"),
 nboot=100, alpha, sigma.sq, theta, theta0, simmap.tree=FALSE, root.age=NULL,
 scaleHeight=FALSE, root.station=TRUE, shift.point=0.5, clade=NULL, mserr="none", 
 diagn=FALSE, quiet=TRUE, warn=FALSE)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data matrix containing species information.}
\item{model}{models to fit to comparative data.}
\item{nboot}{The number of bootstrap replicates.}
\item{alpha}{a numeric vector giving the values of \eqn{\alpha}{alpha} for each selective regime.}
\item{sigma.sq}{a numeric vector giving the values of \eqn{\sigma^2}{sigma^2} for each selective regime.} 
\item{theta}{a numeric vector giving the values of \eqn{\theta}{theta} for each selective regime.}
\item{theta0}{a numeric indicating the starting state, \eqn{\theta_0}{theta_0}}
\item{simmap.tree}{a logical indicating whether the input tree is in SIMMAP format. The default is \code{FALSE}.} 
\item{root.age}{indicates the age of the tree. This is to be used in cases where the "tips" are not contemporary, such as in cases for fossil trees. Default is \code{NULL} meaning latest tip is modern day.}
\item{scaleHeight}{a logical indicating whether the total tree height should be scaled to 1 (see Details). The default is \code{FALSE}.} 
\item{root.station}{a logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated (see Details).} 
\item{shift.point}{the point along a branch where a regime change is assumed to have occurred (if \code{SIMMAP=FALSE}. The default is set to 0.5, or halfway along a branch.}
\item{clade}{a list containing a pair of taxa whose MRCA is the clade of interest.}
\item{mserr}{designates whether a fourth column in the data matrix contains measurement error for each species value ("known"). The measurement error is assumed to be the standard error of the species mean. The default is "none".}
\item{diagn}{a logical indicating whether the full diagnostic analysis should be carried out. The default is \code{FALSE}.}
\item{quiet}{a logical indicating whether progress should be written to the screen. The default is \code{TRUE}.}
\item{warn}{a logical indicating whether a warning should be printed if the number of parameters exceeds ntips/10. The default is \code{FALSE}.}
}
\details{

A simple function for conducting a parametric bootstrap on parameters estimated in OUwie. As before, the input is a tree and a data file. The tree must be of class \dQuote{phylo} and if \code{simmap=FALSE} must contain the ancestral selective regimes as internal node labels. The data file is a dataframe that must have column entries in the following order: [,1] species names and [,2] their current selective regime. The user specifies the simulated parameter values (i.e. \eqn{\alpha}{alpha}, \eqn{\sigma^2}{sigma^2}, \eqn{\theta_0}{theta_0}, \eqn{\theta}{theta}), which is assumed to be the maximum likelihood estimates obtained from an OUwie run.

Note that if \code{root.station} is \code{TRUE} (the default), \eqn{\theta_0}{theta_0} was dropped from the model. In this case, then, \eqn{\theta_0}{theta_0} should be set to the value of the selective regime mapped at the root (i.e., state 1 in the \dQuote{tworegime} example dataset). 

}
\value{

\code{OUwie.boot} returns an object of class \code{OUwie.boot}. This is a matrix of column length equal to the number of parameters, and row length of the number of bootstrap replicates specified.

}
\examples{
\donttest{

data(tworegime)

##First step is estimate parameters under a particular model:
#pp <- OUwie(tree,trait,model=c("OUMV"),root.station=TRUE)

##Second step is to run bootstrap replicates:
#boot.reps <- OUwie.boot(tree,trait,model="OUMV", nboot=10, alpha=pp$solution[1,], 
#sigma.sq=pp$solution[2,],theta=pp$theta[,1], theta0=pp$theta[1,1])

##Finally summarize to obtain the desired confidence -- here is the 95% CI:
#apply(boot.reps, 2, quantile, probs=c(0.025,0.975))

}
}
\references{
Beaulieu J.M., Jhwueng D.C., Boettiger C., and O'Meara B.C. 2012. Modeling stabilizing selection: Expanding the Ornstein-Uhlenbeck model of adaptive evolution. Evolution 66:2369-2383.

O'Meara B.C., Ane C., Sanderson P.C., Wainwright P.C. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60:922-933.

Butler M.A., King A.A. 2004. Phylogenetic comparative analysis: A modeling approach for adaptive evolution. American Naturalist 164:683-695.
}
\author{Jeremy M. Beaulieu}
\keyword{models}
