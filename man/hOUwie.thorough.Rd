\name{hOUwie.thorough}
\alias{hOUwie.thorough}
\title{Rerun a set of hOUwie models with the best mappings of the set.}
\description{This function takes as input a list of hOUwie fits and will automatically model rerun the model set using the best mappings. Given the complexity of a hOUwie model, it is often beneificial to run this to ensure the MLE has been found.}
\usage{
hOUwie.thorough(model.list, ncores = 1)
}
\arguments{
\item{model.list}{A list of model results from several fits of the \code{hOUwie} model. All elements of the list must be of class "houwie".}
\item{ncores}{An integer, specifying the number of cores for running multiple fitting n_starts in parallel. Should generally not exceed the number of CPU cores on a machine, but must be a least 1.}
}
\details{
This function will find the nSim best unique mappings from the input model set and rerun all models based on those mappings alone. Unlike the usual hOUwie, the mappings themselves are not infered and take as given. This makes this function most similar to hOUwie.fixed.
}
\value{
Returns a new model list with maps from the input list.
}
\author{James D. Boyko}
\keyword{hOUwie}
