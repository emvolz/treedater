\name{treedater}
\alias{dater}
\title{
Dating of phylogenetic trees with rate variation.
}
\description{
Estimates the calendar time of the given phylogenetic tree with branches in units of substitutions per site. The calendar time of each sample must also be specified and the length of the sequences used to estimate the tree. 
}
\usage{
dater(tre, sts, s = 1000, omega0 = NA, minblen = NA, maxit = 20, abstol = .01, quiet = FALSE, searchRoot=3)
}

\arguments{
  \item{tre}{
An ape::phylo which describes the phylogeny with branches in units of substitutions per site. Ideally, the tree will be rooted, such as by using an outgroup. If the tree is not rooted, the root position will be estimated using root-to-tip regression. 
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
}

\details{
treedater uses heuristic search to optimise the tmrca's of a phylogeny and the substitution rate. 
An uncorrelated relaxed molecular clock accounts for rate variation between lineages of the phylogeny which is parameterised using a Gamma-Poisson mixture model.
}
\value{
list( omegas = omegas, r = r, theta = gammatheta, Ti = Ti
		 , meanRate = .mean.rate(Ti, r, gammatheta, omegas, td)
		 , loglik = ll )
A list with components
\item{tre }{A phylogenetic tree with branch lengths in calendar time}
\item{meanRate }{The mean substitution rate}
\item{omegas }{The substition rate on each branch of the tree (in same order as tre$edge)}
\item{Ti }{The calendar time of each internal node}
\item{loglik }{The log likelihood of the Gamma-Poisson mixture model of the returned tree}
\item{r }{Gamma shape parameter. If theta = 0, this is the mean rate.}
\item{theta }{Gamma scale parameter}
\item{timeToMRCA}{Estimated time before most recent sample of MRCA}
\item{timeOfMRCA}{Estimated time of MRCA}
}
\references{
Paper in preparation. Please cite github repository.
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
