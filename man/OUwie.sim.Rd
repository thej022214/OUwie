\name{OUwie.sim}
\alias{OUwie.sim}
\title{Generalized Hansen model simulator}
\description{Simulates the Ornstein-Uhlenbeck process of continuous characters evolving under discrete selective regimes.}
\usage{
OUwie.sim(phy=NULL, data=NULL, simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE,
alpha=NULL, sigma.sq=NULL, theta0=NULL, theta=NULL, mserr="none", shift.point=0.5, 
fitted.object=NULL, get.all=FALSE)
}
\arguments{
\item{phy}{a phylogenetic tree, in \code{ape} \dQuote{phylo} format and with internal nodes labeled denoting the ancestral selective regimes}
\item{data}{a dataframe containing species information (see Details). Not necessary to include if \code{simmap=TRUE}.}
\item{simmap.tree}{a logical indicating whether the input tree is in SIMMAP format. The default is \code{FALSE}.}
\item{root.age}{indicates the age of the tree. This is to be used in cases where the "tips" are not contemporary, such as in cases for fossil trees. Default is \code{NULL} meaning latest tip is modern day.}
\item{scaleHeight}{a logical indicating whether the total tree height should be scaled to 1 (see Details). The default is \code{FALSE}.}
\item{alpha}{a numeric vector giving the values of \eqn{\alpha}{alpha} for each selective regime (see Details)}
\item{sigma.sq}{a numeric vector giving the values of \eqn{\sigma^2}{sigma^2} for each selective regime (see Details)}
\item{theta0}{a numeric indicating the starting state, \eqn{\theta_0}{theta_0}}
\item{theta}{a numeric vector giving the values of \eqn{\theta}{theta} for each selective regime (see Details)}
\item{mserr}{designates whether a third column in the data matrix contains measurement error for each species value ("known"). The measurement error is assumed to be the standard error of the species mean. The default is "none".}
\item{shift.point}{the point along a branch where a regime change is assumed to have occurred (if \code{SIMMAP=FALSE}. The default is set to 0.5, or halfway along a branch.}
\item{fitted.object}{a model fit from OUwie to use for simulation.}
\item{get.all}{a logical indicating whether or not the entire simulation history is to be returned. The default is \code{FALSE} meaning that only the tips are returned.}
}
\details{
The input is a tree and a data file OR a fitted OUwie object. The tree must be of class \dQuote{phylo} and if \code{simmap=FALSE} must contain the ancestral selective regimes as internal node labels. The data file is a dataframe that must have column entries in the following order: [,1] species names and [,2] their current selective regime. If \code{mserr="known"} then a third column can be added which contains the measurement error for each species. Note that if \code{simmap=TRUE} no data file is needed. The user specifies the simulated parameter values (i.e. \eqn{\alpha}{alpha}, \eqn{\sigma^2}{sigma^2}, \eqn{\theta_0}{theta_0}, \eqn{\theta}{theta}). Assuming two selective regimes, possible models can be specified as follows (Note that this assumes a stationary distribution at the root):

a. Single rate Brownian motion (\code{BM1}): alpha=c(1e-10,1e-10); sigma.sq=c(0.45,0.45); theta0=1.0; theta=c(0,0).

b. Brownian motion with different rate parameters for each state on a tree (\code{BMS}): alpha=c(1e-10,1e-10) sigma.sq=c(0.45,0.90); theta0=1.0; theta=c(0,0).

c. Ornstein Uhlenbeck with a single optimum for all species (\code{OU1}): alpha=c(0.1,0.1); sigma.sq=c(0.9,0.9); theta0=1; theta=c(1.0,1.0).

d. Ornstein Uhlenbeck model that assumes different state means and a single \eqn{\alpha}{alpha} and \eqn{\sigma^2}{sigma^2} (\code{OUM}): alpha=c(1.0,1.0); sigma.sq=c(0.45,0.45); theta0=1.0; theta=c(1.0,2.0).

e. Ornstein Uhlenbeck model that assumes different state means and multiple \eqn{\sigma^2}{sigma^2} (\code{OUMV}): alpha=c(1.0,1.0); sigma.sq=c(0.45,0.90); theta0=1.0; theta=c(1.0,2.0).

f. Ornstein Uhlenbeck model that assumes different state means and multiple \eqn{\alpha}{alpha} (\code{OUMA}): alpha=c(1.0,0.5); sigma.sq=c(0.45,0.45); theta0=1.0; theta=c(1.0,2.0).

g. Ornstein Uhlenbeck model that assumes different state means and multiple \eqn{\sigma^2}{sigma^2} and \eqn{\alpha}{alpha} (\code{OUMVA}): alpha=c(1.0,0.5); sigma.sq=c(0.45,0.9); theta0=1.0; theta=c(1.0,2.0).

With a fitted OUwie model, it just uses the parameters from that, ignoring any alpha, theta, etc. set in the function.

}
\value{
A dataframe containing, as column entries, [,1] species names, [,2] current regime, [,3] simulated continuous trait, x.
}
\examples{
\donttest{
data(sim.ex)

#Simulate an Ornstein-Uhlenbeck model with different state means
#and a separate alpha and sigma^2 per selective regime
alpha=c(1.0,0.5)
sigma.sq=c(0.45,0.9)
theta0=1.0
theta=c(1.0,2.0)

sim.data<-OUwie.sim(tree,trait,simmap.tree=FALSE,scaleHeight=FALSE,
alpha=alpha,sigma.sq=sigma.sq,theta0=theta0,theta=theta)

#Now fit a model to this and simulate from the fitted results
result <- OUwie(tree, sim.data, model="OUMVA", simmap.tree=FALSE,scaleHeight=FALSE)
sim.data.2 <- OUwie.sim(fitted.object=result)

}
}
\author{Jeremy M. Beaulieu and Brian C. O'Meara}
