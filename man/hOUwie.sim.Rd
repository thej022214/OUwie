\name{hOUwie.sim}
\alias{hOUwie.sim}
\title{Simulate a discrete and continuous character following a Markov and Ornstein-Uhlenbeck model.}
\description{Simulates a discrete and continuous character following the hOUwie model. The function first evolves a discrete character based on the Q matrix provided. Next, it will evolve a continuous character following the given OU parameters and simulated discrete character. Like the hOUwie model, transitions between discrete states are assumed to take place half-way between nodes.}
\usage{
hOUwie.sim(phy, 
           Q, 
           root.freqs, 
           alpha, 
           sigma.sq, 
           theta0, 
           theta)
}
\arguments{
\item{phy}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
\item{Q}{A transition rate matrix with dimensions \code{nStates} by \code{nStates} describing rates of change between discrete states.}
\item{root}{A vector \code{nStates} long with probabilities of the root being in a particular state. For example, for a binary discrete character a root prior of \code{c(1, 0)} would fix the root state in state 1.}
\item{alpha}{A vector \code{nStates} long which gives alpha parameter of the OU model. For a BM model set this to be near 0.}
\item{sigma.sq}{A vector \code{nStates} long which gives the evolutionary rate parameter of the OU model.}
\item{theta0}{A numeric value giving the starting value of the continuous character at the root.}
\item{theta}{A vector \code{nStates} long which gives the phenotypic optima for each regime state.}
}
\details{
The simulation protocol follows the hOUwie model where stochastic maps being produced are based on a node-by-node simulation. I.e., we first simulate the node states given the parameters and then branches are paintined based on transitions occuring half-way between the nodes. 
}
\value{
\item{data}{a dataframe of sp (species), reg (discrete character states), and x (continuous character values).}
\item{simmap}{the history of the discrete character presented as a stochastic map.}
}
\examples{
\donttest{
data(tworegime)
# simulate an OUM model
Q <- matrix(c(-1,1,1,-1), 2, 2)
root.freqs <- c(1, 0)
alpha <- c(2, 2)
sigma.sq <- c(1,1)
theta0 <- 5
theta <- c(5, 10)

simulated_data <- hOUwie.sim(tree, Q, root.freqs, alpha, sigma.sq, theta0, theta)
}
}
\author{James D. Boyko}
\keyword{hOUwie}
