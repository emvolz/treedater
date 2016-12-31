\name{parboot.treedater}
\title{Estimate of confidence intervals using parametric bootstrap for molecular clock dating. }
\description{
	This function simulates phylogenies with branch lengths in units of substitutions per site. Simulations are based on a fitted treedater object which provides parameters of the molecular clock model. The treedater method is applied to each simulated tree providing a Monte Carlo estimate of variance in rates and dates.
}
\usage{
parboot.treedater( td, nreps = 100, level = .95 )
}

\arguments{
\item{td}{
	A fitted treedater object. 
}
\item{nreps}{Integer number of simulations to be carried out}
\item{level}{Confidence region to be computed.}
}

\details{
	If the original treedater fit estimated the root position, root position will also be estimated for each simulation, so the returned trees may have different root positions. Some replicates may converge to a strict clock or a relaxed clock, so the parameter estimates in each replicate may not be directly comparable. It is possible to compute confidence intervals for the times of particular nodes or for estimated sample times by inspecting the output from each fitted treedater object, which is contained in the $trees attribute. 
}

\value{
A list with attributes:
\item{trees}{The fitted treedater objects corresponding to each simulation}
\item{meanRates}{Vector of estimated rates for each simulation}
\item{meanRate_CI}{ Confidence interval for substitution rate}
\item{coef_of_variation_CI}{ Confidence interval for rate variation}
\item{timeOfMRCA_CI}{Confidence interval for time of common ancestor}
}

\references{
Paper in preparation. Please cite github repository.
}
\author{
Erik M Volz
}



\seealso{
dater
}
