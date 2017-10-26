\name{treedater-package}
\alias{treedater-package}
\alias{treedater}
\docType{package}
\title{
\packageTitle{treedater}
}
\description{
\packageDescription{treedater}
}
\details{


\packageDESCRIPTION{treedater}
\packageIndices{treedater}

}
\author{
\packageAuthor{treedater}

Maintainer: \packageMaintainer{treedater}
}
\references{
E.M. Volz and Frost, S.D.W. (2017) Scalable relaxed clock phylogenetic dating. Virus Evolution.
}
\keyword{ package }
\seealso{
}
\examples{

# make a random tree
tre <- rtree(50)
# sample times based on distance from root to tip
sts <- setNames(  dist.nodes( tre)[(length(tre$tip.label)+1), 1:(length(tre$tip.label)+1)], tre$tip.label)
# modify edge length to represent evolutionary distance with rate 1e-3
tre$edge.length <- tre$edge.length * 1e-3
# treedater: 
td <- treedater( tre, sts =sts )


}
