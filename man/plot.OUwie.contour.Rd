\name{plot.OUwie.contour}
\alias{plot.OUwie.contour}
\title{Contour plot}
\description{A plotting function for visualizing likelihood surface for a pair of parameters using OUwie.contour data.}
\usage{
\method{plot}{OUwie.contour}(x, mle.point=NULL, levels=c(0:20*0.1), xlab=NULL, 
ylab=NULL, xlim=NULL, ylim=NULL, col=grey.colors(21, start=0, end=1), ...)
}       
\arguments{
\item{x}{a \code{OUwie.contour} object.}
\item{mle.point}{specifies the color for the maximum likelihood set of parameters. By default the point is not plotted (i.e., set to \code{NULL}.}
\item{levels}{the levels at which to draw contour lines, measured as lnL units away from the best values.}
\item{xlab}{allows users to specify the x-axis label.}
\item{ylab}{allows users to specify the y-axis label.}
\item{xlim}{allows users to specify the lower and upper limits of the x-axis.}
\item{ylim}{allows users to specify the lower and upper limits of the y-axis.}
\item{col}{indicates the color gradient.}
\item{...}{additional parameters to control the plot.}
}
\references{
Beaulieu J.M., Jhwueng D.C., Boettiger C., and O'Meara B.C. 2012. Modeling stabilizing selection: Expanding the Ornstein-Uhlenbeck model of adaptive evolution. Evolution 66:2369-2383.
}
\author{Jeremy M. Beaulieu}
\keyword{plotting}
