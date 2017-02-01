\name{plot.parboot.ltt}
\alias{plot.parboot.ltt}
\title{ Plot estimated lineages through time with CIs based on bootstrap replicates. }
\description{
	This function computes the lineages through time given bootstrap replicate trees. The pseudo-maximum likelihood estimate is plotted alongside CIs based on bootstrap trees. 
}
\usage{
plot.parboot.ltt(pbtd, t0 = NA, res = 100, ... )
}

\arguments{
\item{pbtd}{
	A fitted parboot.treedater object
}
\item{t0}{Lower time limit for LTT plot. If not provided, will use the minimum among bootstrap TMRCAs}
\item{res}{Resolution of time axis}
\item{...}{Additional arguments are passed to ggplot}
}


\value{
A ggplot object. See documentation for ggplot2 package for details of how to show or save these images.
}

\references{
Paper in preparation. Please cite github repository.
}
\author{
Erik M Volz
}



\seealso{
dater, parboot.treedater
}
