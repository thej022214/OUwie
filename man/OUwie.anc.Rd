\name{OUwie.anc}
\alias{OUwie.anc}
\title{Estimate ancestral states given a fitted OUwie model}
\description{Fits ancestral states (a joint estimate) given a fitted OUwie model. Currently, only works for trees with regimes painted on as node labels (rather than a simmap tree). The intended use case is just to visualize what the model is saying about evolution to help intuition (is the model something you can believe in?) rather than a firm estimate you should use to say, Yes, 56.4 MY, the ancestral body size was 17.34 mm. It could likely be anywhere from 1 to 100 mm with about equal chance. Please read details before using this function.}
\usage{
OUwie.anc(fitted.OUwie.object, opts = list("algorithm"="NLOPT_LN_SBPLX",
"maxeval"="1000", "ftol_abs"=0.001), knowledge=FALSE)
}
\arguments{
\item{fitted.OUwie.object}{an object returned from the OUwie function}
\item{opts}{a list of options to pass to nloptr for the ancestral state optimization}
\item{knowledge}{a logical indicating whether or not you have read the documentation The default is FALSE.}
}
\details{
You probably DON'T want to use this function for anything more serious than poking around to make sure the data and model look right (oh, golly: my tips are all between 5 and 15 mm in size, and the ancestral states are 674 mm.). The request to implement ancestral state estimation resulted in many important comments from experts in the public R-SIG-PHYLO discussion forum. Here is a sampling; to see them all, go to \href{https://www.mail-archive.com/r-sig-phylo@r-project.org/msg05257.html}{R-SIG-PHYLO}.

"So in short, yes, you can do it, with any number of methods. But why? If you can answer your biological question with methods that do not involve estimation of a parameter that is inherently fraught with error, it might be better to go another way. Bottom line - use caution and be thoughtful!" -- Marguerite Butler

"I would add an extra caveat to Marguerite's excellent post: Most researchers work with extant taxa only, ignoring extinction. This causes a massive ascertainment bias, and the character states of the extinct taxa can often be very different to the ancestral state reconstructions, particularly if the evolutionary model is wrong. Eg. there has been an evolutionary trend for example. Ancestral state reconstructions based only on extant taxa should be treated as hypotheses to be tested with fossil data. I wouldn't rely on them for much more." -- Simone Blomberg

"While I am at it, let me echo Simone and Marguerite's warnings. The predicted ancestral states will reflect the process you assumed to predict them. Hence, if you  use them to make inferences about evolution, you will recover your own assumptions. I.e. if you predict from a model with no trend, you will find no trend, etc. Many comparative studies are flawed for this reason." -- Thomas Hansen

"Let me add more warnings to Marguerite and Thomas's excellent responses. People may be tempted to infer ancestral states and then treat those inferences as data (and also to infer ancestral environments and then treat those inferences as data). In fact, I wonder whether that is not the main use people make of these inferences. But not only are those inferences very noisy, they are correlated with each other. So if you infer the ancestral state for the clade (Old World Monkeys, Apes) and also the ancestral state for the clade (New World Monkeys, (Old World Monkeys, Apes)) the two will typically not only be error-prone, but will also typically be subject to strongly correlated errors.  Using them as data for further inferences is very dubious. It is better to figure out what your hypothesis is and then test it on the data from the tips of the tree, without the intermediate step of taking ancestral state inferences as observations. The popular science press in particular demands a fly-on-the-wall account of what happened in evolution, and giving them the ancestral state inferences as if they were known precisely is a mistake." -- Joe Felsenstein

"The minor twist I would throw in is that it's difficult to make universal generalizations about the quality of ancestral state estimation.  If one is interested in the ancestral state value at node N, it might be reasonably estimated if it is nested high up within the phylogeny, if the rates of change aren't high, etc. And (local) trends etc might well be reliably inferred.  We are pretty confident that the common ancestor of humans and chimps was larger than many deeper primate ancestors, for instance. If N is the root of your available phylogeny, however, you have to be much more cautious." -- Nick Matzke

"I'll also add that I think there's a great deal to be skeptical of ancestral trait reconstruction even when large amounts of fossil data is available. You can try the exercise yourself: simulate pure BM on a non-ultrametric tree with lots of 'extinct' tips, and you'll still find pretty large confidence intervals on the estimates of the trait values. What does it mean to do ancestral trait reconstruction, if our calculations of uncertainty are that broad?" -- Dave Bapst

These are some of the people who best know the power and limitations of the OU model in phylogenetics. Heed them!

To ensure that you've read this before use, please pass knowledge=TRUE as an argument to the function.
}
\value{
\code{OUwie.anc} returns an object of class \code{OUwie.anc}. This is an OUwie object but with terminal species added at each node representing the ancestral states at each node (which are also included in the data object). There is also a \code{NodeRecon} element in the list that has the optimal ancestral states. There is not currently an estimate of uncertainty, but it is substantial.
}

\examples{
\donttest{
data(tworegime)
fitted <- OUwie(tree,trait,model=c("OUMV"),root.station=TRUE)
recon <- OUwie.anc(fitted, knowledge=FALSE) # read the description above to make this work.
plot(recon)
}
}
\author{Brian C. O'Meara}
