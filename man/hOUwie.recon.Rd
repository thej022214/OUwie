\name{hOUwie.recon}
\alias{hOUwie.recon}
\title{Reconstruct the marginal probability of discrete node states under the hOUwie model.}
\description{Uses the output object from \code{hOUwie} to conduct an ancestral state reconstruction of the discrete regime states.}
\usage{
hOUwie.recon(houwie_obj, 
             nodes      = "all")
}
\arguments{
\item{houwie_obj}{A list of class "houwie".}
\item{nodes}{One of "internal", "external", or "all" to represent reconstructing internal nodes, external (tips) nodes, or all nodes. "External" is only useful if the model has hidden states or there is uncertainty in the tip states.}
}
\details{
Reconstructs the marginal probability of a particular state at a given node. To do this, all possible states at the given node are fixed and the joint likelihoods are evaluated for each possibility. Those joint probabilities are then normalized by the total likelihood so that the probability sums to one.
}
\value{
\item{recon_matrix}{a data.frame containing the marginal probabilities of each state (given by the column).}
}
\references{
Yang Z. 2006. Computational molecular evolution. Oxford University Press Oxford.
}
\author{James D. Boyko}
\keyword{hOUwie}
