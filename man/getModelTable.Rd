\name{getModelTable}
\alias{getModelTable}
\title{Generate a table from a set of hOUwie models describing their relative fit to data.}
\description{This function takes as input a list of hOUwie models and outputs their relative fit to a dataset.}
\usage{
getModelTable(model.list, 
              type = "AIC")
}
\arguments{
\item{model.list}{A list of model results from several fits of the \code{hOUwie} model. All elements of the list must be of class "houwie".}
\item{type}{One of AIC, BIC, or AICc for use during evaluation of relative model fit. AICc is the best option for datasets with a few number of species.}
}
\details{
Models are named based either on the names of the list or the order in which they are provided. If names are given, the rownames will reflect the input. If no names are given, model are assigned names M1 to Mn for the 1 to n models.
}
\value{
\item{model_table}{A data.frame containing the number of parameters (np), total joint log likelihood (lnLik), marginal discrete log likelihood (DiscLik), marginal continuous log likelihood (ContLik), Akaike Information Criterion (AIC), difference in AIC (dAIC), and AIC weight (AICwt).}
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

# get a model table describing the relative fits to the data
getModelTable(model_list)
}
}
}
\references{
Akaike H. 1998. Information Theory and an Extension of the Maximum Likelihood Principle. :15.

Burnham K.P., Anderson D.R. 2002. Model selection and multimodel inference: a practical information-theoretic approach. New York: Springer.
}
\author{James D. Boyko}
\keyword{hOUwie}
