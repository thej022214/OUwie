\name{hOUwie.fixed}
\alias{hOUwie.fixed}
\title{Fit a joint model of discrete and continuous characters via maximum-likelihood with fixed regimes.}
\description{
The hOUwie model jointly estimates the likelihood of a discrete and continuous character by combining the probability of the continuous character given a particular regime and the probability of that discrete regime painting, integrated over many regime paintings. In the standard \code{hOUwie} function, regimes are generated for each iteration of parameters and then evaluated. However, \code{hOUwie.fixed} allows users to provide a map or set of maps which will then be optimized based on their joint probability. Most other inputs and outputs remain the same.

This function takes as main input a list of stochastic maps and a data.frame containing species, character values, and potentially measurement error. The rate category is a way to specify the number of hidden states to be incldued, with one meaning no hidden states present. The structure of the model and association between discrete and continuous characters is specified via discrete_model and continuous model (see Details).
}
\usage{
hOUwie.fixed(simmaps, 
             data, 
             rate.cat, 
             discrete_model, 
             continuous_model, 
             null.model          = FALSE,
             root.p              = "yang", 
             dual                = FALSE, 
             collapse            = TRUE, 
             root.station        = FALSE, 
             get.root.theta      = FALSE, 
             mserr               = "none", 
             lb_discrete_model   = NULL, 
             ub_discrete_model   = NULL, 
             lb_continuous_model = NULL, 
             ub_continuous_model = NULL, 
             recon               = FALSE, 
             nodes               = "internal", 
             p                   = NULL, 
             ip                  = NULL, 
             optimizer           = "nlopt_ln", 
             opts                = NULL, 
             quiet               = FALSE, 
             sample_tips         = FALSE, 
             sample_nodes        = TRUE, 
             adaptive_sampling   = FALSE, 
             diagn_msg           = FALSE,
             make_numeric        = TRUE,
             n_starts            = 1, 
             ncores              = 1)
}

\arguments{
\item{simmaps}{A list of stochastic maps to be evaluated.}
\item{data}{A data frane containing species information. The first column is always species names matching the tip labels of \code{phy}. The \code{2nd to ith} columns are independent discrete characters. The \code{ith+1} column is the continuous character value. Finally, there can be an \code{ith+2} column for measurement error (optional).}
\item{rate.cat}{Specifies the number of rate categories (see Details).}
\item{discrete_model}{Either a user-supplied index of parameters to be optimized or one of "ARD", "SYM", or "ER". ARD: all rates differ. SYM: rates between any two states do not differ. ER: all rates are equal.}
\item{continuous_model}{Either a user-supplied index matrix specifying the continuous model parameters to be estimated or one of "BM1", "BMV", "OU1", "OUA", "OUV", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA" (See also \code{getOUParamStructure}).}
\item{null.model}{A boolean indicating whether the model being run is a character-independent model with rate heterogeneity. Rate.cat must be greater than 1.}
\item{root.p}{a vector used to fix the probabilities at the root, but \dQuote{yang} can also be supplied to use the method of Yang (2006) (see Details).}
\item{dual}{A boolean indicating whether or not to include dual transitions. If true, then transitions two or more states may change at any instant of time. For example, X0Y0 can change directly to X1Y1 without first going through X1Y0 or X0Y1.}
\item{collapse}{A boolean indicating whether to collapse multiple character combinations into only the observed states. For example, if true a two character dataset contained (0,0), (1,0), and (1,1), this would be collapsed into 1,2,3. However, if set to false it would 1,2,4. In combination with a custom rate matrix this allows for the estimation of transitions between unobserved character combinations. Transitions to and from unobserved state combinations may not be estimatable: See Boyko and Beaulieu (2022).} 
\item{root.station}{A boolean indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated (see Details of OUwie function).} 
\item{get.root.theta}{A boolean indicating whether the starting state, \eqn{\theta_0}{theta_0}, should be estimated (see Details of OUwie function).}
\item{mserr}{Designates whether the last column of the data matrix contains measurement error for each species value ("known"). The measurement error is assumed to be the standard error of the species mean. The default is "none".}
\item{lb_discrete_model}{A single value for the lower bound of discrete character evolution during the ML search. The default is \code{1/(Tmax*10000)}.} 
\item{ub_discrete_model}{A single value for the upper bound of discrete character evolution during the ML search. The default is \code{1/(Tmax*0.0001)}.} 
\item{lb_continuous_model}{Three values specifying the lower bound of continuous character evolution during the ML search. Should be provided in the order of lower.bound for alpha, sigma squared, theta. The default is \code{c(1e-10, 1e-10, smallest_continuous_value/10)}.} 
\item{ub_continuous_model}{Three values specifying the upper bound of continuous character evolution during the ML search. Should be provided in the order of upper.bound for alpha, sigma squared, theta. The default is \code{c(log(2)/(0.01 * Tmax), log(2)/(0.01 * Tmax), largest_continuous_value*10)}.} 
\item{recon}{A boolean indicating whether marginal ancestral state reconstruction should take place following the ML search.}
\item{nodes}{If \code{recon=TRUE}, a character or vector indicating which nodes should have their discrete ancestral state characters estimated. If a vector, nodes specified will be reconstructed with 1:Ntip representing tip reconstructions and Ntip+1 representing the root. A user can also indicate "all", "internal", or "external" to reconstruct all nodes, ancestral nodes, or extant nodes respectively. Reconstructing tip states ("external") is only useful if hidden states are included in the model.}
\item{p}{A numeric vector of fixed parameters to be optimized. Parameters are organized such following \code{c(p_trans, p_alpha, p_sigma.squared, p_theta)}}
\item{ip}{A numeric vector indicating initial parameters to start the ML search. Parameters are organized such following \code{c(p_trans, p_alpha, p_sigma.squared, p_theta)}.}
\item{optimizer}{One of "nlopt_ln", "nlopt_gn", or "sann" for a local search of parameter space, global search of paramater space, or simmulated annealing. The default is "nlopt_ln".}
\item{opts}{A list of options to be passed to nloptr.}
\item{quiet}{A boolean indicating whether or not to print messages.}
\item{sample_tips}{A boolean indicating whether or not to sample tip probabilities when generating underlying regimes (see Deatils, but recommended to be set to FALSE).}
\item{sample_nodes}{A boolean indicating whether or not to sample node probabilities when generating underlying regimes (see Deatils, but recommended to be set to TRUE).}
\item{adaptive_sampling}{A boolean indicating whether each iteration will continue sampling stochastic maps until the likelihood of a set of parameters does not improve (recommended TRUE).}
\item{diagn_msg}{A boolean indiciating whether a message of loglikelihoods and parameters should be printed at each iteration. Can be useful for diagnosing problems.}
\item{make_numeric}{A boolean indicating whether the given simmap needs to be converted to numeric map edges matching observed character data. Setting this to FALSE can be useful if you are interested in looking at hidden states since you will can have more mapped states than observed states.}
\item{n_starts}{An integer, specifying the number of trials for fitting the model, using alternative (randomized) starting parameters at each trial. A larger n_starts reduces the risk of landing on a local non-global optimum of the likelihood function, and thus increases the chances of finding the MLE.}
\item{ncores}{An integer, specifying the number of cores for running multiple fitting n_starts in parallel. Should generally not exceed the number of CPU cores on a machine, but must be a least 1.}
}

\details{
The difference between \code{hOUwie.fixed} and \code{hOUwie} is that regimes are provided apriori to hOUwie.fixed and then evaluated. One thing to note is that the discrete probabilities will be treated on a node-by-node basis. That is to say, the exact path probability of the stochastic map is not calculated. Rather, the discrete probability uses the same estimation as hOUwie where we simply evaluate the joint probability of that particular mapping. However, the continuous probability does reflect the exact mapping. This may lead to differences in how hOUwie and hOUwie.fixed preform in empirical situations.

See \code{hOUwie} for general details about the hOUwie model. 
}

\value{
A named list with the following elements:
\item{loglik}{The maximum log-likelihood.}
\item{DiscLik}{The marginal log-likelihood of the discrete state probability for the \code{nSim} maps.}
\item{ContLik}{The marginal log-likelihood of the continuous value probability for the \code{nSim} maps.}
\item{AIC}{Akaike information criterion.}
\item{AICc}{Akaike information criterion corrected for sample-size.}
\item{BIC}{Bayesian information criterion.}
\item{param.count}{Number of parameters in the model.}
\item{solution.disc}{The maximum likelihood estiamte of the transition rate matrix.}
\item{solution.cont}{The maximum likelihood estimate of parameters which describe continuous evolution. Each column corresponds to a discrete state and each row is a different continuous model parameter. Estimated continuous parameters include alpha, sigma.squared, and theta.}
\item{recon}{A data.frame of marginal probabilities of ancestral states given the best fitting model. If ancestral state reconstruction was not specified then \code{NULL}.}
\item{index.disc}{An index matrix specifying which discrete model parameters are estimated.}
\item{index.cont}{An index matrix specifying which continuous model parameters are estimated.}
\item{phy}{User provided phylogenetic tree.}
\item{legend}{A named vector with elements corresponding to input states and names corresponding to a numeric version of those discrete states. Useful for interpreting hOUwie solutions.}
\item{expected_vals}{Expected values given OU parameters averaged over the joint likelihood of the \code{nSim} regime mappings.}
\item{data}{User provided data.frame.}
\item{hOUwie.dat}{Internally used data.frame (provided in case a user is concerned about character matching).}
\item{rate.cat}{The number of rate categories.}
\item{discrete_model}{The user provided index matrix or character string describing the discrete model.}
\item{continuous_model}{The user provided index matrix or character string describing the continuous model.}
\item{root.p}{The root prior used in model estimation.}
\item{time_slice}{Sets how often internodes are sampled. At the moment it is automatically set to be the \code{max(branching.times(phy))+1}.}
\item{root.station}{A boolean indicating whether the starting state, \eqn{\theta_0}{theta_0}, was estimated (recommended FALSE).}
\item{get.root.theta}{A boolean indicating whether the root.theta was included in the model.}
\item{mserr}{Either includes measurement error for each species value ("known"). The default is "none".}
\item{sample_tips}{A boolean indicating whether tips probabilities were sampled each iteration (recommended FALSE).}
\item{sample_nodes}{A boolean indicating whether node probabilities were sampled each iteration (recommended TRUE).}
\item{adaptive_sampling}{A boolean indicating whether each iteration continued sampling stochastic maps until the likelihood of a set of parameters did not improve (recommended TRUE).}
\item{lb_discrete_model}{A vector of lower bounds on discrete character evolution.}
\item{ub_discrete_model}{A vector of upper bounds on discrete character evolution.}
\item{lb_continuous_model}{A vector of lower bounds of alpha, sigma square, and theta.}
\item{ub_continuous_model}{A vector of lower bounds of alpha, sigma square, and theta.}
\item{p}{Vector of maximum likeihood parameters if not given.}
\item{ip}{Vector of starting parameters}
\item{nSim}{Number of stochastic maps for each iteration.}
\item{opts}{List of options input for nloptr optimization.}
\item{quiet}{A boolean indicating whether quiet was TRUE or FALSE.}
\item{all_disc_liks}{A vector of discrete probabilities in the order of the \code{nSim} simmaps.}
\item{all_cont_liks}{A vector of continuous probabilities in the order of the \code{nSim} simmaps.}
\item{simmaps}{\code{nSim} stochastic maps from the best fitting hOUwie parameters (to be removed).}
\item{run_time}{Time to complete the model fit.}
}

\examples{
\donttest{
data(tworegime)

trait[,c(2)] <- c("A", "B")[trait[,c(2)]]
MK <- corHMM:::corHMM(tree, trait[,c(1,2)], 1)
model <- MK$solution
simmap <- corHMM:::makeSimmap(tree=tree, data=trait[,c(1,2)], model=model, 
                              rate.cat=1, nSim=10, nCores=1)

p <- c(0.01664314, 0.01664314, 0.39631218, 0.18684476, 2.2512157, 0.8249509)
pp_cd <- hOUwie.fixed(simmap, trait, rate.cat = 1, discrete_model = "ARD",
                      continuous_model = "OUM", p = p)
# done
}
}

\references{
Beaulieu J.M., Jhwueng D.-C., Boettiger C., O'Meara B.C. 2012. Modeling Stabilizing Selection: Expanding the Ornstein--Uhlenbeck Model of Adaptive Evolution. Evolution. 66:2369--2383.

Beaulieu J.M., O'Meara B.C., Donoghue M.J. 2013. Identifying Hidden Rate Changes in the Evolution of a Binary Morphological Character: The Evolution of Plant Habit in Campanulid Angiosperms. Syst Biol. 62:725--737.

Boyko J.D., Beaulieu J.M. 2021. Generalized hidden Markov models for phylogenetic comparative datasets. Methods Ecol Evol. 12:468--478.

Butler M.A., King A.A. 2004. Phylogenetic Comparative Analysis: A Modeling Approach for Adaptive Evolution. The American Naturalist. 164:683--695.

Hansen T.F. 1997. Stabilizing Selection and the Comparative Analysis of Adaptation. Evolution. 51:1341--1351.

Hansen T.F., Pienaar J., Orzack S.H. 2008. A Comparative Method for Studying Adaptation to a Randomly Evolving Environment. Evolution. 62:1965--1977.

Ho L. si, Ane C. 2014. A Linear-Time Algorithm for Gaussian and Non-Gaussian Trait Evolution Models. Syst Biol. 63:397--408.
}
\author{James D. Boyko}
\keyword{hOUwie}
