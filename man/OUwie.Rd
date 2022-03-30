\name{OUwie}
\alias{OUwie}
\title{Generalized Hansen models}
\description{Fits generalized Ornstein-Uhlenbeck-based Hansen models of continuous characters evolving under discrete selective regimes.}
\usage{
OUwie(phy, data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA",
"TrendyM","TrendyMS"), simmap.tree=FALSE, root.age=NULL,scaleHeight=FALSE,
root.station=FALSE, get.root.theta=FALSE, shift.point=0.5, clade=NULL, mserr="none",
starting.vals=NULL, check.identify=TRUE, algorithm=c("invert", "three.point"),
diagn=FALSE, quiet=FALSE, warn=TRUE, lb = NULL, ub = NULL, opts = list(algorithm = 
"NLOPT_LN_SBPLX", maxeval = "1000", ftol_rel = .Machine$double.eps^0.5))
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data.frame containing species information (see Details).}
\item{model}{models to fit to comparative data (see Details).}
\item{simmap.tree}{a logical indicating whether the input tree is in SIMMAP format. The default is \code{FALSE}.}
\item{root.age}{indicates the age of the tree. This is to be used in cases where the "tips" are not contemporary, such as in cases for fossil trees. Default is \code{NULL} meaning latest tip is modern day.}
\item{scaleHeight}{a logical indicating whether the total tree height should be scaled to 1 (see Details). The default is \code{FALSE}.}
\item{root.station}{a logical indicating whether to assume a random starting point (TRUE) or a fixed starting point (FALSE) (see Details).}
\item{get.root.theta}{a logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated (see Details).}
\item{shift.point}{the point along a branch where a regime change is assumed to have occurred (if \code{SIMMAP=FALSE}. The default is set to 0.5, or halfway along a branch.}
\item{clade}{a list containing a pair of taxa whose MRCA is the clade of interest (see Details).}
\item{mserr}{designates whether a fourth column in the data matrix contains measurement error for each species value ("known"). The measurement error is assumed to be the standard error of the species mean. The default is "none".}
\item{starting.vals}{a vector of initial values for the optimization search. For OU models, two must be supplied, with the first being the initial alpha value and the second being the initial sigma squared. For BM models, just a single value is needed.}
\item{check.identify}{a logical indicating whether to check that the user-supplied regime paintings will produce identifiable theta estimates. The default is \code{TRUE}.}
\item{algorithm}{designates whether the standard matrix inversion ('invert') or the faster 'three-point' algorithm of Ho and Ane (2013) should be used.}
\item{diagn}{a logical indicating whether the full diagnostic analysis should be carried out. The default is \code{FALSE}.}
\item{quiet}{a logical indicating whether progress should be written to the screen. The default is \code{FALSE}.}
\item{warn}{a logical indicating whether a warning should be printed if the number of parameters exceeds ntips/10. The default is \code{TRUE}.}
\item{lb}{if algorithm == "invert" a single value indicating the lower bound for the parameter values. if algorithm == "three.point", a vector of length 3 with position 1 as alpha lower bound, position 2 as sigma.sq's lower bound, position 3 as theta's lower bound. when set to NULL it will be a default value of 1e-9. Note that even if the model you're using doesn't include alpha (e.g. BM1), it must be included in this vector, but it will not be used to set any bounds.}
\item{ub}{if algorithm == "invert" a single value indicating the upper bound for the parameter values. if algorithm == "three.point", a vector of length 3 with position 1 as alpha upper bound, position 2 as sigma.sq's upper bound, position 3 as theta's upper bound. when set to NULL it will be a default value of 100. Note that even if the model you're using doesn't include alpha (e.g. BM1), it must be included in this vector, but it will not be used to set any bounds.}
\item{opts}{a list of options to pass to nloptr for the optimization: useful to adjust for faster, coarser searches}
}
\details{
This function fits various likelihood models for continuous characters evolving under discrete selective regimes. The function returns parameter estimates and their approximate standard errors. The R package \code{nloptr} provides a common interface to NLopt, an open-source library for nonlinear optimization. The likelihood function is maximized using the bounded subplex optimization routine (\code{NLOPT_LN_SBPLX}). As input all \code{OUwie} requires is a tree and a trait data.frame. The tree must be of class \dQuote{phylo} and must contain the ancestral selective regimes as internal node labels. Internal node labels can be applied manually or from some sort of ancestral state reconstruction procedure (BayesTraits, \code{ape}, \code{diversitree}, SIMMAP, etc.), which would then be brought into OUwie. This is essentially what is required by \code{ouch} and Brownie (though Brownie provides built-in ancestral state reconstruction capabilities). The trait data.frame must have column entries in the following order: [,1] species names, [,2] current selective regime, and [,3] the continuous trait of interest. Alternatively, if the user wants to incorporate measurement error (\code{mserr}="known"), then a fourth column, [,4] must be included that provides the standard error estimates for each species mean. However, a global measurement error for all taxa can be estimated from the data (\code{mserr}="est"); is not well tested, so use at your own risk. Also, a user can specify a particular clade as being in a different selective regime, by specifying a pair of species whose mrca is the root of the clade of interest [e.g., \code{clade}=c("taxaA","taxaB")]. OUwie will automatically assign internal node labels and update the data matrix according to this clade designation.

Possible models are as follows: single-rate Brownian motion (\code{model=BM1}), Brownian motion with different rate parameters for each state on a tree (\code{model=BMS}), Ornstein-Uhlenbeck model with a single optimum for all species (\code{model=OU1}), Ornstein-Uhlenbeck model with different state means and a single \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} acting all selective regimes (\code{model=OUM}), and new Ornstein-Uhlenbeck models that assume different state means as well as either multiple \eqn{\sigma^2}{sigma^2} (\code{model=OUMV}), multiple \eqn{\alpha}{alpha} (\code{model=OUMA}), or multiple \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} per selective regime (\code{model=OUMVA}).

By default, we drop the root optima and absorb the weight into whatever regime the root is in. In previous version we used to incorrectly refer to this as "stationarity". True stationarity assumes that the starting state comes from a distribution, and the covariance requires an additional variance term to account for the fact that, up until T=0, the lineage is assumed to have been evolving in the ancestral regime. We have added this in for the OU1 and OUM models only (root.station=TRUE).

Note, too, that when specifying the BMS model also be mindful of the root.station flag. When root.station=FALSE, the non-censored model of O'Meara et al. 2006 is invoked (i.e., a single regime at the root is estimated), and when root.station==TRUE the group mean model of Thomas et al. 2006 (i.e., the number of means equals the number of regimes). The latter case appears to be a strange special case of OU, in that it behaves similarly to the OUMV model, but without selection. I would say that this is more consistent with the censored test of O'Meara et al. (2006), as opposed to having any real connection to OU. In any case, more work is clearly needed to understand the behavior of the group means model, and therefore, I recommend setting root.station=FALSE in the BMS case.

The Hessian matrix is used as a means to estimate the approximate standard errors of the model parameters and to assess whether they are the maximum likelihood estimates. The variance-covariance matrix of the estimated values of \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} are computed as the inverse of the Hessian matrix and the standard errors are the square roots of the diagonals of this matrix. The Hessian is a matrix of second-order derivatives and is approximated in the R package \code{numDeriv}. So, if changes in the value of a parameter results in sharp changes in the slope around the maximum of the log-likelihood function, the second-order derivative will be large, the standard error will be small, and the parameter estimate is considered stable. On the other hand, if the second-order derivative is nearly zero, then the change in the slope around the maximum is also nearly zero, indicating that the parameter value can be moved in any direction without greatly affecting the log-likelihood. In such situations, the standard error of the parameter will be large.

For models that allow \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} to vary (i.e., \code{OUMV}, \code{OUMA}, and \code{OUMVA}), the complexity of the model can often times be greater than the information that is contained within the data. As a result one or many parameters are poorly estimated, which can cause the function to return a log-likelihood that is suboptimal. This has great potential for poor model choice and incorrect biological interpretations. An eigendecomposition of the Hessian can provide an indication of whether the search returned the maximum likelihood estimates. If all the eigenvalues of the Hessian are positive, then the Hessian is positive definite, and all parameter estimates are considered reliable. However, if there are both positive and negative eigenvalues, then the objective function is at a saddlepoint and one or several parameters cannot be estimated adequately. One solution is to just fit a simpler model. Another is to actually identify the offending parameters. This can be done through the examination of the eigenvectors. The row order corresponds to the entries in \code{index.matrix}, the columns correspond to the order of values in \code{eigval}, and the larger the value of the row entry the greater the association between the corresponding parameter and the eigenvalue. Thus, the largest values in the columns associated with negative eigenvalues are the parameters that are causing the objective function to be at a saddlepoint.
}
\value{
\code{OUwie} returns an object of class \code{OUwie}. This is a list with elements:
\item{$loglik}{the maximum log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$BIC}{Bayesian information criterion.}
\item{$mBIC}{modified Bayesian information criterion of Ho and Ane (2014).}
\item{$model}{The model being fit.}
\item{$param.count}{The number of parameters counted in the model.}
\item{$solution}{a matrix containing the maximum likelihood estimates of \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2}.}
\item{$theta}{a matrix containing the maximum likelihood estimates of \eqn{\theta}{theta} and its standard error.}
\item{$solution.se}{a matrix containing the approximate standard errors of \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2}. The standard error is calculated as the diagonal of the inverse of the Hessian matrix.}
\item{$tot.state}{A vector of names for the different regimes}
\item{$index.mat}{The indices of the parameters being estimated are returned. The numbers correspond to the row in the \code{eigvect} and can useful for identifying the parameters that are causing the objective function to be at a saddlepoint (see Details).}
\item{$simmap.tree}{A logical indicating whether the input phylogeny is a SIMMAP formatted tree.}
\item{$root.age}{The user-supplied age at the root of the tree.}
\item{$split.point}{The user-supplied point at which regime changes are assumed to have occurred.}
\item{$opts}{Internal settings of the likelihood search.}
\item{$data}{User-supplied dataset.}
\item{$phy}{User-supplied tree.}
\item{$root.station}{A logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, was estimated.}
\item{$starting.vals}{A vector of user-supplied initial search parameters.}
\item{$lb}{The lower bound set.}
\item{$ub}{The upper bound set.}
\item{$iterations}{Number of iterations of the likelihood search that were executed.}
\item{$get.root.theta}{Indicates whether the root.theta was included in the model.}
\item{$regime.weights}{A table containing parameter estimates and the weights for time spent in each regime for each tip.}
\item{$eigval}{The eigenvalues from the decomposition of the Hessian of the likelihood function. If any \code{eigval<0} then one or more parameters were not optimized during the likelihood search (see Details).}
\item{$eigvect}{The eigenvectors from the decomposition of the Hessian of the likelihood function is returned (see Details).}
\item{$new.start}{The vector of values to use if you want to restart the run from this point (starting.vals for a new run).}
\item{$algorithm}{The algorithm used to estimate parameters}
}
\examples{
\donttest{
data(tworegime)

#Plot the tree and the internal nodes to highlight the selective regimes:
select.reg<-character(length(tree$node.label))
select.reg[tree$node.label == 1] <- "black"
select.reg[tree$node.label == 2] <- "red"
plot(tree)
nodelabels(pch=21, bg=select.reg)

}

\dontrun{
#To see the first 5 lines of the data matrix to see what how to
#structure the data:
trait[1:5,]

#Now fit an OU model that allows different sigma^2:
OUwie(tree,trait,model=c("OUMV"))

#Fit an OU model based on a clade of interest:
OUwie(tree,trait,model=c("OUMV"), clade=c("t50", "t64"), algorithm="three.point")

#For large trees, it may be useful to have ways to restart the search (due to
#finite time per run on a computing cluster, for example). You can do this
#by changing settings of OUwie runs. For example:

run1 <- OUwie(tree,trait,model=c("OUMV"), root.station=FALSE, algorithm="invert", 
opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="500", "ftol_abs"=0.001))

save(run1, file="run1.rda")

#Then, later or in a different session:

load("run1.rda")

run2 <- OUwie(tree, trait, model=c("OUMV"), algorithm="three.point", 
opts = list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="500", "ftol_abs"=0.001),
starting.vals=run1$new.start)

#run2 will start off where run1 stopped.
}
}
\references{
Beaulieu J.M., Jhwueng D.C., Boettiger C., and O'Meara B.C. 2012. Modeling stabilizing selection: Expanding the Ornstein-Uhlenbeck model of adaptive evolution. Evolution 66:2369-2383.

O'Meara B.C., Ane C., Sanderson P.C., Wainwright P.C. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60:922-933.

Butler M.A., King A.A. 2004. Phylogenetic comparative analysis: A modeling approach for adaptive evolution. American Naturalist 164:683-695.

Ho, L.S.T., and C. Ane. 2014. Intrinsic inference difficulties for trait evolution with Ornstein-Uhlenbeck models. Methods in Ecology and Evolution, 5: 1133-1146.

Thomas G.H., Freckleton R.P., and Szekely T. 2006. Comparative analysis of the influence of developmental mode on phenotypic diversification rates in shorebirds. Proceedings of the Royal Society, B. 273:1619-1624.
}
\author{Jeremy M. Beaulieu and Brian C. O'Meara}
\keyword{models}
