\name{OUwie.joint}
\alias{OUwie.joint}
\title{Joint optimization for multiple traits}
\description{Fits generalized Ornstein-Uhlenbeck-based Hansen models for multiple continuous characters evolving under discrete selective regimes.}
\usage{
OUwie.joint(phy, data, model=c("BM1","BMS","OU1","OUM","OUMV",
"OUMVr","OUMA","OUMAr","OUMVA","OUMVAr"), ntraits, allfree=TRUE, 
simmap.tree=FALSE, scaleHeight=FALSE, root.station=TRUE, 
lb=0.000001, ub=1000, mserr="none", diagn=FALSE, quiet=FALSE)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data matrix containing species information (see Details).}
\item{model}{models to fit to comparative data (see Details).}
\item{ntraits}{a numeric indicating the number of traits to analyze.}
\item{allfree}{a logical indicating whether the parameters vary across the multiple traits (see Details). The default is \code{TRUE}.}
\item{simmap.tree}{a logical indicating whether the input tree is in SIMMAP format. The default is \code{FALSE}.} 
\item{scaleHeight}{a logical indicating whether the total tree height should be scaled to 1 (see Details). The default is \code{FALSE}.} 
\item{root.station}{a logical indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated (see Details).} 
\item{lb}{lower bound for the likelihood search. The default is \code{lb=0.000001}. Note that this value must be greater than zero.}
\item{ub}{upper bound for the likelihood search. The default is \code{ub=1000}.}
\item{mserr}{designates whether the data matrix contains measurement error for each species value ("known"). The default is "none".}
\item{diagn}{a logical indicating whether the full diagnostic analysis should be carried out. The default is \code{FALSE}.}
\item{quiet}{a logical indicating whether progress should be written to the screen. The default is \code{FALSE}.}
}
\details{
This function jointly estimates parameters under various likelihood models for several continuous characters evolving under discrete selective regimes. As with the regular \code{OUwie} function all that is required is a tree and a data matrix. The tree must be of class \dQuote{phylo} and must contain the ancestral selective regimes as internal node labels. Internal node labels can be applied manually or from some sort of ancestral state reconstruction procedure (BayesTraits, \code{ape}, \code{diversitree}, SIMMAP, etc.), which would then be brought into OUwie. The data matrix must have column entries in the following order: [,1] species names, [,2] current selective regime, with all other columns containing the set of continuous traits to optimize jointly. Currently, this function does not allow for measurement error to be incorporated, but future versions will (\code{mserr}="known"). Finally, a user can simply specify a particular clade as being in a different selective regime, by inputting a pair of species whose mrca is the root of the clade of interest [e.g., \code{clade}=c("taxaA","taxaB")]. OUwie will automatically assign internal node labels and update the data matrix according to this clade designation. 

Possible models are the same as before: single-rate Brownian motion (\code{model=BM1}), Brownian motion with different rate parameters for each state on a tree (\code{model=BMS}), Ornstein-Uhlenbeck model with a single optimum for all species (\code{model=OU1}), Ornstein-Uhlenbeck model with different state means and a single \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} acting all selective regimes (\code{model=OUM}), and new Ornstein-Uhlenbeck models that assume different state means as well as either multiple \eqn{\sigma^2}{sigma^2} (\code{model=OUMV}), multiple \eqn{\alpha}{alpha} (\code{model=OUMA}), or multiple \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} per selective regime (\code{model=OUMVA}). We also allow for several constrained models: a model that assumes different state means as well as multiple \eqn{\sigma^2}{sigma^2} but keeps \eqn{\alpha}{alpha} constant across all traits (\code{model=OUMVr}); a model that assumes different state means as well as either multiple \eqn{\alpha}{alpha} but keeps \eqn{\sigma^2}{sigma^2} constant across all traits (\code{model=OUMAr}). Finally, when allfree=FALSE, the model specified will be considered global and applied to all traits rather than estimating the same model separately for each trait. 
}
\value{
\code{OUwie.joint} returns an object of class \code{OUwie.joint}. This is a list with elements:
\item{$loglik}{the maximum negative log-likelihood.}
\item{$AIC}{Akaike information criterion.}
\item{$AICc}{Akaike information criterion corrected for sample-size.}
\item{$model}{The model being fit}
\item{$solution}{a matrix containing the maximum likelihood estimates of \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2}.}
\item{$theta}{a matrix containing the maximum likelihood estimates of \eqn{\theta}{theta} and it standard error.}
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
\item{$ntraits}{Number of traits optimized}
}
\examples{
#data(tworegime)
#pp <- OUwie(tree,trait,model=c("OUMV"),root.station=TRUE)

#trait1 <- OUwie.sim(tree,trait, alpha=pp$solution[1,], 
#sigma.sq=pp$solution[2,],theta=pp$theta[,1], theta0=pp$theta[1,1])
#trait2 <- OUwie.sim(tree,trait, alpha=pp$solution[1,], 
#sigma.sq=pp$solution[2,],theta=pp$theta[,1], theta0=pp$theta[1,1])
#trait.multi <- cbind(trait[,1:2],trait1[,3],trait2[,3])

##Fit a global OUMV model applied to all traits:
#oumv.global <- OUwie.joint(tree,trait.multi, model="OUMV", ntraits=2, allfree=FALSE)

##Now fit a OUMV model applied to each separately, but optimized jointly:
#oumv.joint <- OUwie.joint(tree,trait.multi, model="OUMV", ntraits=2, allfree=TRUE)
}
\references{

Leslie, A.B., Beaulieu, J.M., Crane, P.R., and Donoghue, M.J. 2014. Cone size is related to branching architecture in conifers. New Phytologist doi:10.111/nph.12864.

Beaulieu J.M., Jhwueng D.C., Boettiger C., and O'Meara B.C. 2012. Modeling stabilizing selection: Expanding the Ornstein-Uhlenbeck model of adaptive evolution. Evolution 66:2369-2383.

}
\author{Jeremy M. Beaulieu}
\keyword{models}
