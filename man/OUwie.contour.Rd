\name{OUwie.contour}
\alias{OUwie.contour}
\title{Generates data for contour plot of likelihood surface}
\description{Generates the likelihood surface for pairs of free parameters for generalized Ornstein-Uhlenbeck-based Hansen models of continuous characters evolving under discrete selective regimes.}
\usage{
OUwie.contour(OUwie.obj, focal.params=c("alpha_1", "sigma.sq_1"), 
focal.params.lower=c(0,0), focal.params.upper=c(5,5), nreps=1000, n.cores=NULL)
}       
\arguments{
\item{OUwie.obj}{an object of class \dQuote{OUwie} that contains the focal parameters for conducting the likelihood surface search.}
\item{focal.params}{a vector specifying the parameters that you would like the likelihood surface. The format is parameter_regime -- that is, for theta in regime 1, the input would be "theta_1".} 
\item{focal.params.lower}{a vector specifying the lower bounds for the parameters. The values need to be in the order as the focal parameters.}
\item{focal.params.upper}{a vector specifying the upper bounds for the parameters. The values need to be in the order as the focal parameters.}
\item{nreps}{the number points to use to estimate the likelihood surface (see Details).}
\item{n.cores}{specifies the number of independent processors to conduct the analysis. The default is \code{NULL}.}
}
\details{
This function samples a set of points to estimate the likelihood surface for any pair of parameters, letting the other parameters find their own optima. This process can be very slow, as it involves optimization \code{nrep} times (though with two fewer parameters than with the chosen model, as the focal parameter values are fixed). It uses a latin hypercube design to sample points across the user-defined range of the focal parameters.

The pair of parameters to examine is passed by focal.param. The parameters need to be one of three: theta, alpha, sigma.sq. For example, to do a plot of sigma.sq from the first regime and alpha from the second regime, one would pass \code{focal.param = c( "sigma.sq_1", "alpha_2")}. As another example, if the regimes are characters like, flower color, the focal parameter would be \code{focal.param = c( "sigma.sq_Red", "sigma.sq_Blue")}. 

This returns a data.frame with the last two columns being the values of the points examined and the first column the loglik of those points. The first row contains the MLE. The data.frame can be incorporated into a plotting function to obtain a contour plot (see plot.OUwie.contour).
}
\value{
\item{surface.data}{the parameter values and loglik}
\item{focal.params}{the vector specifying the parameter pair for which likelihood surface is evaluated}
\item{focal.params.lower}{the vector specifying the lower bounds for the parameter pair.}
\item{focal.params.upper}{the vector specifying the uuper bounds for the parameter pair.}
}
\references{
Beaulieu J.M., Jhwueng D.C., Boettiger C., and O'Meara B.C. 2012. Modeling stabilizing selection: Expanding the Ornstein-Uhlenbeck model of adaptive evolution. Evolution 66:2369-2383.
}
\author{Jeremy M. Beaulieu}

