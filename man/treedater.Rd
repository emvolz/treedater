\name{treedater}
\alias{dater}
\title{
Dating of phylogenetic trees with rate variation.
}
\description{
Estimates the calendar time of nodes in the given phylogenetic tree with branches in units of substitutions per site. The calendar time of each sample must also be specified and the length of the sequences used to estimate the tree. If the tree is not rooted, this function will estimate the root position. 
}
\usage{
dater(tre, sts, s = 1000, omega0 = NA, minblen = NA, maxit = 100, abstol = .001, quiet = TRUE, searchRoot=5, temporalConstraints=TRUE, strictClock=FALSE, estimateSampleTimes=NULL, estimateSampleTimes_densities = list())
}

\arguments{
  \item{tre}{
An ape::phylo which describes the phylogeny with branches in units of substitutions per site. This may be a rooted or unrooted tree. If unrooted, the root position will be estimated by checking multiple candidates chosen by root-to-tip regression.  If the tree has multifurcations, these will be resolved and a binary tree will be returned.
}
  \item{sts}{
Vector of sample times for each tip in phylogenetic tree. Vector must be named with names corresponding to tre$tip.label. 
}
  \item{s}{
Sequence length. 
}
 \item{omega0}{
Initial guess of the mean substitution rate (substitutions per site per unit time). If not provided, will guess using root to tip regression. 
}
 \item{minblen}{
Minimum branch length in calendar time. By default, this will be the range of sample times (max - min) divided by sample size. 
}
 
  \item{maxit}{
Maximum number of iterations.
}

\item{abstol}{
	Difference in log likelihood between successive iterations for convergence. 
}

\item{searchRoot}{
	Will search for the optimal root position using the top matches from root-to-tip regression. 
	If searchRoot=x, dates will be estimated for x trees, and the estimate with the highest likelihood will be returned.
}

\item{quiet}{
	If TRUE, will suppress messages during execution.
}

\item{temporalConstraints}{
	If TRUE, will enforce the condition that an ancestor node in the phylogeny occurs before all progeny. Equivalently, this will preclude negative branch lengths. Note that execution is faster if this option is FALSE. 
}

\item{strictClock}{
	If TRUE, will fit a Poisson evolutionary model without rate variation. 
}

\item{estimateSampleTimes}{
	If some sample times are not known with certainty, bounds can be provided with this option. This should take the form of a data frame with columns 'lower' and 'upper' providing the sample time bounds for each uncertain tip. Row names of the data frame should correspond to elements in tip.label of the input tree. Tips with sample time bounds in this data frame do not need to appear in the *sts* argument, however if they are included in *sts*, that value will be used as a starting condition for optimisation. 
}

\item{estimateSampleTimes_densities}{
	An optional named list of log densities which would be used as priors for unknown sample times. Names should correspond to elements in tip.label with uncertain sample times. 
}

}

\details{
treedater uses heuristic search to optimise the tmrca's of a phylogeny and the substitution rate. 
An uncorrelated relaxed molecular clock accounts for rate variation between lineages of the phylogeny which is parameterised using a Gamma-Poisson mixture model.
}
\value{
A phylo object which is rooted and with time-scaled branch lengths. 
There are additional data elements: 
\item{meanRate }{The estimated mean substitution rate}
\item{omegas }{The substition rate on each branch of the tree (in same order as tre$edge)}
\item{Ti }{The calendar time of each internal node}
\item{loglik }{The log likelihood of the Gamma-Poisson mixture model of the returned tree}
\item{r }{Gamma shape parameter. If theta = 0, this is the mean rate.}
\item{theta }{Gamma scale parameter}
\item{timeToMRCA}{Estimated time before most recent sample of MRCA}
\item{timeOfMRCA}{Estimated time of MRCA}
\item{coef_of_variation} Coefficient of variation of estimated rates
\item{clock} 'strict' or 'relaxed'
}
\references{
E.M. Volz and Frost, S.D.W. (2017) Scalable relaxed clock phylogenetic dating. Virus Evolution.
}
\author{
Erik M Volz
}



\seealso{
ape::chronos
}
\examples{

# make a random tree
tre <- rtree(50)
# sample times based on distance from root to tip
sts <- setNames(  dist.nodes( tre)[(length(tre$tip.label)+1), 1:(length(tre$tip.label)+1)], tre$tip.label)
# modify edge length to represent evolutionary distance with rate 1e-3
tre$edge.length <- tre$edge.length * 1e-3
# treedater: 
td <- dater( tre, sts =sts )

}
\keyword{ misc }
