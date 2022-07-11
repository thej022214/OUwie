\name{getOUParamStructure}
\alias{getOUParamStructure}
\title{Generate a continuous model parameter structure}
\description{Given a specific model type (see Details) and number of observed states create a matrix which describes the continuous model structure for use in \code{hOUwie}. This function can be used to generate a set of character-dependent and character-independent models for multi-model inferences.}
\usage{
getOUParamStructure(model, 
                    nObsState,
                    rate.cat   = 1, 
                    null.model = FALSE)
}
\arguments{
\item{model}{One of "BM1", "BMV", "OU1", "OUA", "OUV", "OUM", "OUVA", "OUMV", "OUMA", "OUMVA".}
\item{nObsState}{An integer specifying the total number of observed states.}
\item{rate.cat}{An integer specifying the number of rate categories/ hidden states. If \code{null.model = TRUE}, then this should be at least 2 to allow for character-independent rate heterogeneity.}
\item{null.model}{A boolean indicating whether the model is character-independent with rate heterogeneity (\code{TRUE}) or character-dependent (\code{FALSE}).}
}
\details{
The basic models structures supported and their primary reference (so far as I know) can be found below. Additionally, many intermediate models can be formed by users who modify the default rate matricies when there are multiple rate categories or discrete states. E.g., it may be that states 1 and 3 share an optimum that is unique from state 2. This model is not automatically created by \code{getOUParamStructure}, but can be easily created by the user.

BM1: Brownian motion model with a single evolutionary rate (Felsenstein 1985).

BMV: Brownian motion model with multiple, regime dependent, evolutionary rates (O'Meara et al. 2006).

OU1: Ornstein-Uhlenbeck model with a single optima, alpha, and sigma.sq (Hansen 1997).

OUA: Ornstein-Uhlenbeck model with a single optima, variable alpha, and single sigma.sq (Boyko  et al. 2022).

OUV: Ornstein-Uhlenbeck model with a single optima, single alpha, and variable sigma.sq (Boyko  et al. 2022).

OUM: Ornstein-Uhlenbeck model with variable optima, single alpha, and single sigma.sq (Butler and King 2004).

OUVA: Ornstein-Uhlenbeck model with single optima, variable alpha, and variable sigma.sq (Boyko et al. 2022).

OUMV: Ornstein-Uhlenbeck model with variable optima, single alpha, and variable sigma.sq (Beaulieu et al. 2012).

OUMA: Ornstein-Uhlenbeck model with variable optima, variable alpha, and single sigma.sq (Beaulieu et al. 2012).

OUMVA: Ornstein-Uhlenbeck model with variable optima, variable alpha, and variable sigma.sq (Beaulieu et al. 2012).
}
\value{
Returns a matrix which describes the continuous model structure for use in \code{hOUwie}. 
}
\examples{
\donttest{
\dontrun{
# for the following example we will imagine we have a single discrete character with 3 states.
# these matrices could be assigned to a variable and input into the continuous_model arguement
# in hOUwie for later model comparison.

# there are several model structures but they can be broken into:
# Character dependent (CD) which 
getOUParamStructure("OUMA", 3, 1, FALSE)
# Character independent (CID) which can have rate heterogeneity 
getOUParamStructure("OUMA", 3, 2, TRUE)
# Hybrid model (HYB) which has both character dependent and character independent rate heterogeneity
getOUParamStructure("OUMA", 3, 2, FALSE)

# from the different uses of the function above, notice how the parameters are associated 
# with particular regimes to create the different classes of model. CD has a separate 
# parameter for each state, whereas CID always has the different observed states associated 
# with the same parameter value. Obviously, different mixes of these are allowed and the 
# hybrid model is the most general form which can then be further constrained.
}
}
}
\references{
Beaulieu J.M., Jhwueng D.-C., Boettiger C., O'Meara B.C. 2012. Modeling Stabilizing Selection: Expanding the Ornstein--Uhlenbeck Model of Adaptive Evolution. Evolution. 66:2369--2383.

Boyko J.D., O'Meara B.C., Beaulieu J.M. 2022. In prep.

Butler M.A., King A.A. 2004. Phylogenetic Comparative Analysis: A Modeling Approach for Adaptive Evolution. The American Naturalist. 164:683--695.

Felsenstein J. 1985. Phylogenies and the Comparative Method. Am. Nat. 125:1--15.

Hansen T.F. 1997. Stabilizing Selection and the Comparative Analysis of Adaptation. Evolution. 51:1341--1351.

O'Meara B.C., Ane C., Sanderson M.J., Wainwright P.C. 2006. Testing for Different Rates of Continuous Trait Evolution Using Likelihood. Evolution. 60:922--933.
}
\author{James D. Boyko}
\keyword{model structure}
