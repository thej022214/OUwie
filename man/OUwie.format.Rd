\name{OUwie.format}
\alias{OUwie.format}
\title{Format data and tree for OUwie}
\description{Simplifies conversion of a tree and, optionally, data, for OUwie.}
\usage{
OUwie.format(phy, tip.regimes=NULL, tip.data=NULL, tip.fog=NULL, 
tip.fog.percentage=NULL, verbose=TRUE)
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and optionally with internal nodes labeled denoting the ancestral selective regimes}
\item{tip.regimes}{a named vector or one-column data.frame with rownames with the regimes for each tip}
\item{tip.data}{a named vector or one-column data.frame with rownames with the trait values for each tip}
\item{tip.fog}{a named vector or one-column data.frame with rownames with the tip fog for each tip}
\item{tip.fog.percentage}{a value for the percentage of the trait value to use as the tip fog}
\item{verbose}{a logical indicating whether or not to print messages}
}
\details{
OUwie expects a tree with regimes mapped on nodes or a simmap tree. It also wants a data.frame with a column for taxon names, one for regime, one (JUST one) for trait value, and optionally one for tip fog. This function gets all of those formatted.

If you pass in no data, just the tree, the function will return a tree with regime 1 appearing everywhere, as well as a trait data.frame with regime 1 and trait value NA everywhere. It also includes zero tip fog for all taxa. This is very unrealistic, of course, and creates a bias in results (any variance in tips must be caused by the variance in the evolutionary process if tip fog is assumed to be zero). 

For `tip.regimes`, `tip.data`, and `tip.fog` the user can pass in a vector with the values and the names corresponding to the names of the tips on the tree or a single column data.frame where the rownames match the taxon names on the tree. The names of the input data and the names on the tree do not need to be in the same order, but they do need to match: "Homo sapiens" and "Homo_sapiens" do not match, for example.

Another way of handling tip dfog is assuming it's a fixed percentage of the tip value: "My uncertainty in squid arm lengths is 67 percent".  To do this, pass in the percentage: `tip.fog.percentage=67`. Note that it is per_cent: that is, for 67 percent put in 67, not 0.67. It will only work if there is already a measurement for the trait.
}
\value{
A list with the tree ('tree') and the data ('data')
}

\author{Brian C. O'Meara}
