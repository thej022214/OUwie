\name{fix.kappa}
\alias{fix.kappa}
\title{Adjust tree for matrix condition}
\description{Iteratively deletes taxa with shortest tip length to try to get a variance covariance matrix with good matrix condition.}
\usage{
fix.kappa(phy, data, threshold = log(40))
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data.frame containing species information (see Details).}
\item{threshold}{log(condition), as measured by kappa(), which is too large}
}
\details{
Internally, OUwie uses an algorithm that can perform poorly when the variance covariance matrix is poorly conditioned (which can happen if two columns are very similar, as when the divergence depth of two species is very recent). This does not mean there is anything wrong with the biology, just that the numerical algorithms perform poorly in that case. If it's a model that can be fit in phylolm or geiger, those packages use a differnt algorithm that is more robust to this. What this function does is take your original tree and data and deletes taxa with the shortest branches, in order, to try to get a starting tree with generally good condition. Deleting data is always a sad thing, but this can result in a more accurate estimate of the likelihood and parameter values.
}
\value{
This returns a list with two elements:
\item{$phy}{the phylogeny with taxa deleted.}
\item{$data}{the data with taxa deleted.}
}
\author{Brian C. O'Meara}
\keyword{kappa}
