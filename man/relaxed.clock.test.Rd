\name{relaxed.clock.test}
\alias{relaxed.clock.test}
\title{Use parametric bootstrap to test if relaxed clock offers improved fit to data.}
\description{
	This function simulates phylogenies with branch lengths in units of substitutions per site. Simulations are based on a fitted treedater object which provides parameters of the molecular clock model. The coefficient of variation of rates is estimated using a relaxed clock model applied to strict clock simulations. Estimates of the CV is then compared to the null distribution provided by simulations. 
}
\usage{
relaxed.clock.test( ... )
}

\arguments{
\item{...}{
	The arguments are the same as for the dater function. 
}
}

\details{
	This function will print the optimal clock model under the test and the distribution of rate CV under the null hypothesis (strict clock). 
}

\value{
A dater object under either a strict or relaxed clock. 
}

\references{
E.M. Volz and Frost, S.D.W. (2017) Scalable relaxed clock phylogenetic dating. Virus Evolution.
}
\author{
Erik M Volz & Simon DW Frost
}



\seealso{
dater
}
