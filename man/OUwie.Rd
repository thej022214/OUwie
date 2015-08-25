\name{OUwie}
\alias{OUwie}
\title{Generalized Hansen models}
\description{Fits generalized Ornstein-Uhlenbeck-based Hansen models of continuous characters evolving under discrete selective regimes.}
\usage{
OUwie(phy, data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA",
"TrendyM","TrendyMS"), simmap.tree=FALSE, scaleHeight=FALSE, root.station=TRUE,
clade=NULL, mserr="none", diagn=FALSE, quiet=FALSE, warn=TRUE)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data matrix containing species information (see Details).}
\item{model}{models to fit to comparative data (see Details).}
\item{simmap.tree}{a logical indicating whether the input tree is in SIMMAP format. The default is \code{FALSE}.} 
\item{scaleHeight}{a logical indicating whether the total tree height should be scaled to 1 (see Details). The default is \code{FALSE}.} 
\item{root.station}{a logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated (see Details).} 
\item{clade}{a list containing a pair of taxa whose MRCA is the clade of interest (see Details).}
\item{mserr}{designates whether a fourth column in the data matrix contains measurement error for each species value ("known"). The measurement error is assumed to be the standard error of the species mean. The default is "none".}
\item{diagn}{a logical indicating whether the full diagnostic analysis should be carried out. The default is \code{FALSE}.}
\item{quiet}{a logical indicating whether progress should be written to the screen. The default is \code{FALSE}.}
\item{warn}{a logical indicating whether a warning should be printed if the number of parameters exceeds ntips/10. The default is \code{TRUE}.}
}
\details{
This function fits various likelihood models for continuous characters evolving under discrete selective regimes. The function returns parameter estimates and their approximate standard errors. The R package \code{nloptr} provides a common interface to NLopt, an open-source library for nonlinear optimization. The likelihood function is maximized using the bounded subplex optimization routine (\code{NLOPT_LN_SBPLX}). As input all \code{OUwie} requires is a tree and a data matrix. The tree must be of class \dQuote{phylo} and must contain the ancestral selective regimes as internal node labels. Internal node labels can be applied manually or from some sort of ancestral state reconstruction procedure (BayesTraits, \code{ape}, \code{diversitree}, SIMMAP, etc.), which would then be brought into OUwie. This is essentially what is required by \code{ouch} and Brownie (though Brownie provides built-in ancestral state reconstruction capabilities). The data matrix must have column entries in the following order: [,1] species names, [,2] current selective regime, and [,3] the continuous trait of interest. Alternatively, if the user wants to incorporate measurement error (\code{mserr}="known"), then a fourth column, [,4] must be included that provides the standard error estimates for each species mean. However, a global measurement error for all taxa can be estimated from the data (\code{mserr}="est"); is not well tested, so use at your own risk. Also, a user can specify a particular clade as being in a different selective regime, by inputting a pair of species whose mrca is the root of the clade of interest [e.g., \code{clade}=c("taxaA","taxaB")]. OUwie will automatically assign internal node labels and update the data matrix according to this clade designation. 

The initial implementation followed \code{ouch} in that the tree is automatically rescaled so that the branch lengths were in proportion to the total height of the tree. However, this makes the results inconsistent with other implementations such as Brownie or \code{geiger}. Therefore, we have allow the user to choose whether the tree should be rescaled or not. Note that the when \code{scaleHeight=FALSE} the bounds will have to be adjusted to the appropriate scale.

Possible models are as follows: single-rate Brownian motion (\code{model=BM1}), Brownian motion with different rate parameters for each state on a tree (\code{model=BMS}), Ornstein-Uhlenbeck model with a single optimum for all species (\code{model=OU1}), Ornstein-Uhlenbeck model with different state means and a single \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} acting all selective regimes (\code{model=OUM}), and new Ornstein-Uhlenbeck models that assume different state means as well as either multiple \eqn{\sigma^2}{sigma^2} (\code{model=OUMV}), multiple \eqn{\alpha}{alpha} (\code{model=OUMA}), or multiple \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} per selective regime (\code{model=OUMVA}).

If \code{root.station} is \code{TRUE} (the default), \eqn{\theta_0}{theta_0} is dropped from the model. Under these conditions it is assumed that the starting value is distributed according to the stationary distribution of the OU process. This would not fit a biological scenario involving moving away from an ancestral state, but it does fit a scenario of evolution at a steady state. Dropping \eqn{\theta_0}{theta_0} from the model can sometimes stabilize estimates of the primary optima, especially in situations where the estimates of \eqn{\theta}{theta} in the full model are non-sensical. In regards to the accuracy of estimating \eqn{\theta_0}{theta_0}, it is important to note that in simulation, as \eqn{\alpha}{alpha} increases estimates of \eqn{\theta_0}{theta_0} converge to zero. Thus, when \eqn{\alpha}{alpha} is large (i.e. \eqn{\alpha}{alpha}>2) it is likely that any inference of an evolutionary trend will be an artifact and positively misleading.

The Hessian matrix is used as a means to estimate the approximate standard errors of the model parameters and to assess whether they are the maximum likelihood estimates. The variance-covariance matrix of the estimated values of \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} are computed as the inverse of the Hessian matrix and the standard errors are the square roots of the diagonals of this matrix. The Hessian is a matrix of second-order derivatives and is approximated in the R package \code{numDeriv}. So, if changes in the value of a parameter results in sharp changes in the slope around the maximum of the log-likelihood function, the second-order derivative will be large, the standard error will be small, and the parameter estimate is considered stable. On the other hand, if the second-order derivative is nearly zero, then the change in the slope around the maximum is also nearly zero, indicating that the parameter value can be moved in any direction without greatly affecting the log-likelihood. In such situations, the standard error of the parameter will be large. 

For models that allow \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} to vary (i.e., \code{OUMV}, \code{OUMA}, and \code{OUMVA}), the complexity of the model can often times be greater than the information that is contained within the data. As a result one or many parameters are poorly estimated, which can cause the function to return a log-likelihood that is suboptimal. This has great potential for poor model choice and incorrect biological interpretations. An eigendecomposition of the Hessian can provide an indication of whether the search returned the maximum likelihood estimates. If all the eigenvalues of the Hessian are positive, then the Hessian is positive definite, and all parameter estimates are considered reliable. However, if there are both positive and negative eigenvalues, then the objective function is at a saddlepoint and one or several parameters cannot be estimated adequately. One solution is to just fit a simpler model. Another is to actually identify the offending parameters. This can be done through the examination of the eigenvectors. The row order corresponds to the entries in \code{index.matrix}, the columns correspond to the order of values in \code{eigval}, and the larger the value of the row entry the greater the association between the corresponding parameter and the eigenvalue. Thus, the largest values in the columns associated with negative eigenvalues are the parameters that are causing the objective function to be at a saddlepoint.
}
\value{
\code{OUwie} returns an object of class \code{OUwie}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$model}{The model being fit}
\item{$solution}{a matrix containing the maximum likelihood estimates of \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2}.}
\item{$theta}{a matrix containing the maximum likelihood estimates of \eqn{\theta}{theta} and it standard error.}
\item{$solution.se}{a matrix containing the approximate standard errors of \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2}. The standard error is calculated as the diagonal of the inverse of the Hessian matrix.}
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

#Plot the tree and the internal nodes to highlight the selective regimes:
select.reg<-character(length(tree$node.label)) 
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"
plot(tree) 
nodelabels(pch=21, bg=select.reg)

#Not run
#To see the first 5 lines of the data matrix to see what how to
#structure the data:
#trait[1:5,]

#Now fit an OU model that allows different sigma^2: 
#OUwie(tree,trait,model=c("OUMV"),root.station=TRUE)

#Fit an OU model based on a clade of interest:
#OUwie(tree,trait,model=c("OUMV"), root.station=TRUE, clade=c("t50", "t64"))

}
\references{
Beaulieu J.M., Jhwueng D.C., Boettiger C., and O'Meara B.C. 2012. Modeling stabilizing selection: Expanding the Ornstein-Uhlenbeck model of adaptive evolution. Evolution 66:2369-2383.

O'Meara B.C., Ane C., Sanderson P.C., Wainwright P.C. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60:922-933.

Butler M.A., King A.A. 2004. Phylogenetic comparative analysis: A modeling approach for adaptive evolution. American Naturalist 164:683-695.
}
\author{Jeremy M. Beaulieu and Brian C. O'Meara}
\keyword{models}
