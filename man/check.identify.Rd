\name{check.identify}
\alias{check.identify}
\title{A test of regime identifiability}
\description{Ho and Ane test for determining whether all regimes form connected subtrees, making both the ancestral state and the regime optima unidentifiable.}
\usage{
check.identify(phy, data, simmap.tree=FALSE, quiet=FALSE)
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes.}
\item{data}{a data.frame containing species information (see Details).}
\item{simmap.tree}{a logical indicating whether the input tree is in SIMMAP format. The default is \code{FALSE}.}
\item{quiet}{a logical indicating whether messages should be written to the screen. The default is \code{FALSE}.}
}
\value{
This returns a vector with two elements, with the first being an indicator of identifiability (0=unidentifiable, 1=identifiable, and if \code{get.penalty=TRUE}, the second is the penalty used for the modified BIC.}
\references{
Ho, L.S.T., and C. Ane. 2014. Intrinsic inference difficulties for trait evolution with Ornstein-Uhlenbeck models. Methods in Ecology and Evolution, 5: 1133-1146.
}
\author{Jeremy M. Beaulieu and Brian C. O'Meara}
