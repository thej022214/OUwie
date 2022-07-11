\name{hOUwie}
\alias{hOUwie}
\title{Fit a joint model of discrete and continuous characters via maximum-likelihood.}
\description{
The hOUwie model jointly estimates the likelihood of a discrete and continuous character by combining the probability of the continuous character given a particular regime and the probability of that discrete regime painting, integrated over many regime paintings. Specifically, we combine hidden Markov models of discrete character evolution (Beaulieu et al. 2013; Boyko and Beaulieu 2021) with generalized Ornstein-Uhlenbeck models (Hansen 1997; Butler and King 2004; Hansen et al. 2008; Beaulieu et al. 2012; Ho and Ane 2014a). 

This function takes as main input a phylogenetic tree and a data.frame containing species, character values, and potentially measurement error. The rate category is a way to specify the number of hidden states to be incldued, with one meaning no hidden states present. The structure of the model and association between discrete and continuous characters is specified via discrete_model and continuous model (see Details).}
\usage{
hOUwie(phy, 
       data, 
       rate.cat, 
       discrete_model, 
       continuous_model, 
       null.model          = FALSE,
       nSim                = 100, 
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
       adaptive_sampling   = TRUE, 
       diagn_msg           = FALSE,
       n_starts            = 1, 
       ncores              = 1)
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{A data frane containing species information. The first column is always species names matching the tip labels of \code{phy}. The \code{2nd to ith} columns are independent discrete characters. The \code{ith+1} column is the continuous character value. Finally, there can be an \code{ith+2} column for measurement error (optional).}
\item{rate.cat}{Specifies the number of rate categories (see Details).}
\item{discrete_model}{Either a user-supplied index of parameters to be optimized or one of "ARD", "SYM", or "ER". ARD: all rates differ. SYM: rates between any two states do not differ. ER: all rates are equal.}
\item{continuous_model}{Either a user-supplied index matrix specifying the continuous model parameters to be estimated or one of "BM1", "BMV", "OU1", "OUA", "OUV", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA" (See also \code{getOUParamStructure}).}
\item{null.model}{A boolean indicating whether the model being run is a character-independent model with rate heterogeneity. Rate.cat must be greater than 1.}
\item{nSim}{The number of stochastic maps evaluated per iteration of the ML search.} 
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
\item{ip}{A numeric vector indicating initial parameters to start the ML search. Parameters are organized such following \code{c(p_trans, p_alpha, p_sigma.squared, p_theta)}}
\item{optimizer}{One of "nlopt_ln", "nlopt_gn", or "sann" for a local search of parameter space, global search of paramater space, or simmulated annealing. The default is "nlopt_ln".}
\item{opts}{A list of options to be passed to nloptr.}
\item{quiet}{A boolean indicating whether or not to print messages.}
\item{sample_tips}{A boolean indicating whether or not to sample tip probabilities when generating underlying regimes (see Deatils, but recommended to be set to FALSE).}
\item{sample_nodes}{A boolean indicating whether or not to sample node probabilities when generating underlying regimes (see Deatils, but recommended to be set to TRUE).}
\item{adaptive_sampling}{A boolean indicating whether each iteration will continue sampling stochastic maps until the likelihood of a set of parameters does not improve (recommended TRUE).}
\item{diagn_msg}{A boolean indiciating whether a message of loglikelihoods and parameters should be printed at each iteration. Can be useful for diagnosing problems.}
\item{n_starts}{An integer, specifying the number of trials for fitting the model, using alternative (randomized) starting parameters at each trial. A larger n_starts reduces the risk of landing on a local non-global optimum of the likelihood function, and thus increases the chances of finding the MLE.}
\item{ncores}{An integer, specifying the number of cores for running multiple fitting n_starts in parallel. Should generally not exceed the number of CPU cores on a machine, but must be a least 1.}
}
\details{
Our model is composed of two processes: one that describes the evolution of a discrete character and the other describes the evolution of a continuous character. To model the evolution of a single continuous character we use an Ornstein-Uhlenbeck (OU) model (Hansen 1997; Butler and King 2004; Hansen et al. 2008; Beaulieu et al. 2012; Ho and Ane 2014a). This model combines the stochastic evolution of a trait through time with a deterministic component which models the tendency for a trait to evolve towards an optimum. In this model, the value of a trait, \code{X_t}, is pulled towards an optimum, \code{theta}, at a rate scaled by the parameter \code{alpha}. The optimum, \code{theta}, is a piecewise constant on intervals and takes values in a finite set \code{theta}. This can represent the set of "selective regimes", "regimes", or Simpson's "adaptive zones" (Cressler et al. 2015), though it is consistent with a variety of true underlying microevolutionary models (Hansen 2014). Additionally, random deviations are introduced by Gaussian white noise dB(t), which is distributed as a normal random variable with mean zero and variance equal to \code{sigma^2}. Thus, \code{sigma^2} is a constant describing the rate of stochastic evolution around the optimum. We use the set of extensions introduced by Beaulieu et al. (2012). This allows for multiple primary optima \code{theta} in which both the pull strength (\code{alpha}) and the rate of stochastic evolution (\code{sigma^2}) can vary across the phylogeny. 

Our model allows for multiple rate categories (specified by the rate.cat argument). The default is one rate class in which only observed discrete states evolve. However, there are two main reasons one may be interested in increasing the number of rate categories. First, rate heterogeneity throughout the tree of life is more of a rule than a possibility. Not all things will evolve in the same way at the same time. For example, woody and herbaceous plants may evolve one way on the mainland and a completely different way on island. The challenge then is to allow for heterogeneity in the evolutionary process when we do not know what the 'mainland' or 'island' variables are. This is what hidden Markov models allow us to do (Boyko and Beaulieu 2021). The second reason to include hidden rate categories is to reduce the error rates of finding a false correlation. This has been discussed elsewhere in depth (see Maddison and FitzJohn 2015, Uyeda et al. 2018, Boyko and Beaulieu 2022). The problem lies in if we compare models with rate heterogeneity to models without it. Most character dependent models (correlation models) allow for different ways for the characters to evolve. In fact, their reliable inference depends on being able to define the differences between, for example, body size evolution on islands and mainlands. However, most character-independent models have no way to allow for variable ways for body size to evolve. So, whether or not the evolution of body size is connected to island systems, we may be biased towards selecting correlation models simply because they allow body size to have different rates of evolution. 

It is not uncommon to map a discrete character on a phylogeny using a model of evolution for that discrete character and then run various continuous trait models using that mapping (or a stochastic set of these mappings), and OUwie has ways to use such mappings. However, that mapping is completely uninformed by the continuous trait -- it could be that a nearly identical mapping fits the continuous trait far better. A different approach is to do the regime detection using the continuous trait only, as in SURFACE or OUwie.dredge. But if there is a reasonable discrete character (perhaps with hidden states), it can make the most sense to get the likelihood integrating across all different discrete mappings, summing the likelihood for the discrete and continuous traits while doing so. Because looking at all possible reconstructions is infeasible, and looking at them from a distribution from stochastic mapping biases it towards discrete-only mappings, we approximate this by summing over a large number of potential mappings centered on those that seem a good fit to the data.

One way that hOUwie differs from other implementations of automatic regime detection is that we explicitly calculate the the probability of discrete characters (D) and stochastic mapping (z). The difficulties of calculating the pathway probabilities explicitly are outlined in Boyko et al. (2022), but ultimately led us to use an approximation. Our approximation relies on a finite number of degree-2 internodes and uses the standard Chapman-Kolmgorov equation to calculate the probabilities of beginning in a particular state i and ending in state j (Pagel 1994) and is identical to a joint probability of a set of state reconstructions (Yang 2006). 

\code{discrete_model} defines the model between discrete states. This option controls the number of independent rates and the correspondence between independent and dependent rates. If a character, then it must be one of "ER", "SYM", "ARD". If an integer matrix, then it defines a custom parametric structure for the transition rates, by mapping entries of the transition matrix to a set of independent transition-rate parameters (numbered 1,2, and so on). All cells with the number 3 represent transitions that have the same estimated rate. Additionally, ambiguities (polymorphic taxa or taxa missing data) are assigned likelihoods following Felsenstein (2004, p. 255). Taxa with missing data are coded \dQuote{?} with all states observed at a tip. Polymorphic taxa are coded with states separated by an \dQuote{&}. For example, if a trait has four states and taxon A is observed to be in state 1 and 3, the character would be coded as \dQuote{1&3}. Missing data are treated as ambiguous for all states, thus all states for taxa missing data are assigned a likelihood of 1.0 (some software may incorrectly assign a likelihood of 1.0/number of possible states; note this if comparing likelihoods). For example, for a four-state character (i.e. DNA), a taxon missing data will have likelihoods of all four states equal to 1.0 [e.g. L(A)=1.0, L(C)=1.0, L(G)=1.0, L(T)=1.0].

\code{continuous_model} option controls the number of independent continuous model parameters and the correspondence between continuous model parameters and underlying regime. If a character, then it must be one of "BM1", "BMV", "OU1", "OUM", "OUA", "OUV", "OUMV", "OUMA", "OUVA", or "OUMVA". These correspond to allowing theta ("OUM"), sigma_square ("BMV", "OUV"), or alpha ("OUA") to vary. If an integer matrix, then it defines a custom parametric structure for the transition rates, by mapping entries of the continuous model parameters to the states of the underlying regimes. Note that when some taxa have phenotypic values below 0, we add 50 to all continuous trait values to assist the optimization, but the estimated parameter values will be scaled back to the original values.

For each set of parameters evaluated during the maximum likelihood search, a set of stochastic mappings are generated to evaluate the discrete and continuous likelihoods. To do this efficiently, we first approximate of the conditional state probabilities at nodes. This is what happens if \code{sample_nodes=TRUE}. The conditional state probability, unlike the more common marginal reconstruction or joint state reconstruction, calculates the probability that a node has a particular state value conditioned only on the observations of its descendants. For a particular focal node, we calculate the probability of the observing all pairwise descendant values given the OU model parameters, integrated over all possible rootward node states, and observed tipward discrete states. This process is repeated for each iteration of parameters if \code{adaptive_sample=TRUE}, in which case it is repeated until the joint likelihood for a given set of mappings does not improve.

The algorithm used to calculate the likelihood described in Beaulieu et al. (2012) involves matrix inversion -- a computationally costly procedure. Therefore, we implement a linear-time computation of the likelihood of Gaussian trait models following (Ho and Ane 2014). 

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
require(OUwie)
data(tworegime)

##First step is to define our discrete model:
# we use corHMM functions to help create rate matricies
ER_1 <- corHMM:::getStateMat4Dat(data = trait[,c(1,2)], model = "ER")$rate.mat
# if we wanted to make an independent model we could use more corHMM functions
rate_class_matrix <- corHMM:::getRateCatMat(2)
# this gives us a two rate class matrix
ER_2 <- corHMM:::equateStateMatPars(corHMM:::getFullMat(list(ER_1, ER_1), 
                                                        rate_class_matrix), c(1:4))

##Second step is to define our continuous model:
# character dependent
OUM_CD <- getOUParamStructure(model = "OUM", nObsState = 2, rate.cat = 1, null.model = FALSE)
# character independent
OUM_CID <- getOUParamStructure(model = "OUM", nObsState = 2, rate.cat = 2, null.model = TRUE)


##Third step is to fit our model (here we are using fixed params from a previous ML search):
p <- c(0.01664314, 0.39631218, 0.18684476, 2.25121568, 0.82495093) # MLE
pp_cd <- hOUwie(tree, trait, rate.cat = 1, discrete_model = ER_1, 
                continuous_model = OUM_CD, nSim = 25, p = p)
p <- c(0.008956173, 0.013314238, 0.142922735, 43.563960366, 1.308719702) # MLE
pp_cid <- hOUwie(tree, trait, rate.cat = 2, discrete_model = ER_2, 
                 continuous_model = OUM_CID, nSim = 25, p = p)

## We can compare this to the result from OUwie on a fixed map to see a difference 
# in estimates of theta
pp2 <- OUwie(tree,trait,model=c("OUM"),root.station=FALSE, algorithm="three.point")

pp_cd$solution.cont
pp2$solution
}
}
\references{
Beaulieu J.M., Jhwueng D.-C., Boettiger C., O'Meara B.C. 2012. Modeling Stabilizing Selection: Expanding the Ornstein--Uhlenbeck Model of Adaptive Evolution. Evolution. 66:2369--2383.

Beaulieu J.M., O'Meara B.C., Donoghue M.J. 2013. Identifying Hidden Rate Changes in the Evolution of a Binary Morphological Character: The Evolution of Plant Habit in Campanulid Angiosperms. Syst Biol. 62:725--737.

Boyko J.D., Beaulieu J.M. 2021. Generalized hidden Markov models for phylogenetic comparative datasets. Methods Ecol Evol. 12:468--478.

Butler M.A., King A.A. 2004. Phylogenetic Comparative Analysis: A Modeling Approach for Adaptive Evolution. The American Naturalist. 164:683--695.

Cressler C.E., Butler M.A., King A.A. 2015. Detecting Adaptive Evolution in Phylogenetic Comparative Analysis Using the Ornstein--Uhlenbeck Model. Systematic Biology. 64:953--968.

Felsenstein J. 2004. Inferring phylogenies. Sunderland MA: Sinauer Associates.

Hansen T.F. 1997. Stabilizing Selection and the Comparative Analysis of Adaptation. Evolution. 51:1341--1351.

Hansen T.F., Pienaar J., Orzack S.H. 2008. A Comparative Method for Studying Adaptation to a Randomly Evolving Environment. Evolution. 62:1965--1977.

Hansen T.F. 2014. Use and Misuse of Comparative Methods in the Study of Adaptation. In: Garamszegi L.Z., editor. Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology. Berlin, Heidelberg: Springer Berlin Heidelberg. p. 351--379.

Ho L. si, Ane C. 2014. A Linear-Time Algorithm for Gaussian and Non-Gaussian Trait Evolution Models. Syst Biol. 63:397--408.
}
\author{James D. Boyko}
\keyword{hOUwie}
