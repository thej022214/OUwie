\name{getModelAvgParams}
\alias{getModelAvgParams}
\title{Model average the parameter estimates over severl hOUwie fits.}
\description{This function takes as input a list of hOUwie fits and will automatically model average the parameter estimates. It returns model averaged parameter values.}
\usage{
getModelAvgParams(model.list, 
                  BM_alpha_treatment = "zero", 
                  type               = "BIC",
                  force              = TRUE)
}
\arguments{
\item{model.list}{A list of model results from several fits of the \code{hOUwie} model. All elements of the list must be of class "houwie".}
\item{BM_alpha_treatment}{A boolean indicating whether alpha from BM models should be treated as zero or those models should be removed when averaging alpha values. Currently zero is the only allowed value.}
\item{type}{One of AIC, BIC, or AICc for use during evaluation of relative model fit. AICc is the best option for datasets with a few number of species.}
\item{force}{A boolean indicating whether to force potentially failed model fits to be included in the model averaging.}
}
\details{
The AIC weight of each model is first evaluated and used as the weighting for each models' contribution to the overall parameter value (Burnham and Anderson 2002). Tip values are estimated based on the weighted average of the joint probability of the fitted stochastic maps.
}
\value{
Returns a named list with the following elements:
\item{weighted_tip_values}{Model averaged parameter estimate for each tip.}
}
\examples{
\donttest{
\dontrun{
# fit several possible models (we're using fixed parameters here)
p <- c(0.01664314, 0.39631218, 0.18684476, 2.25121568, 0.82495093) # MLE
Model_X <- hOUwie(tree, trait, rate.cat = 1, discrete_model = "ER", 
                  continuous_model = "OUM", nSim = 25, p = p)
p <- c(0.01664314, 0.39631218, 0.18684476, 0.25, 2.25121568, 0.82495093) # not MLE
Model_Y <- hOUwie(tree, trait, rate.cat = 1, discrete_model = "ER", 
                  continuous_model = "OUMV", nSim = 25, p = p)
p <- c(0.01664314, 0.39631218, 0.09631218, 0.18684476, 2.25121568, 0.82495093) # not MLE
Model_Z <- hOUwie(tree, trait, rate.cat = 1, discrete_model = "ER", 
                  continuous_model = "OUMA", nSim = 25, p = p)

# put the model results into a list
model_list <- list(Model_X = Model_X, Model_Y = Model_Y, Model_Z = Model_Z)

# model average the parameters
getModelAvgParams(model_list)
}
}
}
\references{
Beaulieu J.M., Jhwueng D.-C., Boettiger C., O'Meara B.C. 2012. Modeling Stabilizing Selection: Expanding the Ornstein--Uhlenbeck Model of Adaptive Evolution. Evolution. 66:2369--2383.

Burnham K.P., Anderson D.R. 2002. Model selection and multimodel inference: a practical information-theoretic approach. New York: Springer.

Butler M.A., King A.A. 2004. Phylogenetic Comparative Analysis: A Modeling Approach for Adaptive Evolution. The American Naturalist. 164:683--695.

Hansen T.F. 1997. Stabilizing Selection and the Comparative Analysis of Adaptation. Evolution. 51:1341--1351.
}
\author{James D. Boyko}
\keyword{hOUwie}
