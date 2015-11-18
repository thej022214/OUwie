\name{OUwie.sim}
\alias{OUwie.sim}
\title{Generalized Hansen model simulator}
\description{Simulates the Ornstein-Uhlenbeck process of continuous characters evolving under discrete selective regimes.}
\usage{
OUwie.sim(phy, data=NULL, simmap.tree=FALSE, scaleHeight=FALSE, alpha, 
sigma.sq, theta0, theta, mserr="none")
}       
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes}
\item{data}{a dataframe containing species information (see Details). Not necessary to include if \code{simmap=TRUE}.}
\item{simmap.tree}{a logical indicating whether the input tree is in SIMMAP format. The default is \code{FALSE}.} 
\item{scaleHeight}{a logical indicating whether the total tree height should be scaled to 1 (see Details). The default is \code{FALSE}.}
\item{alpha}{a numeric vector giving the values of \eqn{\alpha}{alpha} for each selective regime (see Details)}
\item{sigma.sq}{a numeric vector giving the values of \eqn{\sigma^2}{sigma^2} for each selective regime (see Details)} 
\item{theta0}{a numeric indicating the starting state, \eqn{\theta_0}{theta_0}}
\item{theta}{a numeric vector giving the values of \eqn{\theta}{theta} for each selective regime (see Details)}
\item{mserr}{designates whether a third column in the data matrix contains measurement error for each species value ("known"). The measurement error is assumed to be the standard error of the species mean. The default is "none".}
}
\details{
The input is a tree and a data file. The tree must be of class \dQuote{phylo} and if \code{simmap=FALSE} must contain the ancestral selective regimes as internal node labels. The data file is a dataframe that must have column entries in the following order: [,1] species names and [,2] their current selective regime. If \code{mserr="known"} then a third column can be added which contains the measurement error for each species. Note that if \code{simmap=TRUE} no data file is needed. The user specifies the simulated parameter values (i.e. \eqn{\alpha}{alpha}, \eqn{\sigma^2}{sigma^2}, \eqn{\theta_0}{theta_0}, \eqn{\theta}{theta}). Assuming two selective regimes, possible models can be specified as follows (Note that this assumes a stationary distribution at the root):

a. Single rate Brownian motion (\code{BM1}): alpha=c(1e-10,1e-10); sigma.sq=c(0.45,0.45); theta0=1.0; theta=c(0,0).

b. Brownian motion with different rate parameters for each state on a tree (\code{BMS}): alpha=c(1e-10,1e-10) sigma.sq=c(0.45,0.90); theta0=1.0; theta=c(0,0).

c. Ornstein Uhlenbeck with a single optimum for all species (\code{OU1}): alpha=c(0.1,0.1); sigma.sq=c(0.9,0.9); theta0=1; theta=c(1.0,1.0).

d. Ornstein Uhlenbeck model that assumes different state means and a single \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} (\code{OUM}): alpha=c(1.0,1.0); sigma.sq=c(0.45,0.45); theta0=1.0; theta=c(1.0,2.0).

e. Ornstein Uhlenbeck model that assumes different state means and multiple \eqn{\sigma^2}{sigma^2} (\code{OUMV}): alpha=c(1.0,1.0); sigma.sq=c(0.45,0.90); theta0=1.0; theta=c(1.0,2.0).

f. Ornstein Uhlenbeck model that assumes different state means and multiple \eqn{\alpha}{alpha} (\code{OUMA}): alpha=c(1.0,0.5); sigma.sq=c(0.45,0.45); theta0=1.0; theta=c(1.0,2.0).

g. Ornstein Uhlenbeck model that assumes different state means and multiple \eqn{\sigma^2}{sigma^2} and \eqn{\alpha}{alpha} (\code{OUMVA}): alpha=c(1.0,0.5); sigma.sq=c(0.45,0.9); theta0=1.0; theta=c(1.0,2.0).

}
\value{
A dataframe containing, as column entries, [,1] species names, [,2] current regime, [,3] simulated continuous trait, x.
}
\examples{
data(sim.ex)

#Simulate an Ornstein-Uhlenbeck model with different state means
#and a separate alpha and sigma^2 per selective regime
alpha=c(1.0,0.5)
sigma.sq=c(0.45,0.9)
theta0=1.0
theta=c(1.0,2.0)

sim.data<-OUwie.sim(tree,trait,simmap.tree=FALSE,scaleHeight=FALSE,
alpha=alpha,sigma.sq=sigma.sq,theta0=theta0,theta=theta)
}
\author{Jeremy M. Beaulieu and Brian C. O'Meara}
