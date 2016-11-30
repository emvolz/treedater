#treedater
`treedater` estimates the calendar time of the given phylogenetic tree with branches in units of substitutions per site. The calendar time of each sample must also be specified and the length of the sequences used to estimate the tree. 

`treedater` uses heuristic search to optimise the tmrca's of a phylogeny and the substitution rate. 
An uncorrelated relaxed molecular clock accounts for rate variation between lineages of the phylogeny which is parameterised using a Gamma-Poisson mixture model.

To cite: A paper is in preparation. Please cite github repository:
* Volz, E.M. (2016). treedater [software].


## basic usage
```{r}
dater( tre, sts, s)
```
where 
* `tre` is an `ape::phylo` phylogeny, 
* `sts` is a named vector of sample times for each tip in `tre`
* `s` is the length of the genetic sequences used to estimate `tre`

## simple example
```{r}
require(treedater)
# make a random tree
tre <- rtree(50)
# sample times based on distance from root to tip
sts <- setNames(  dist.nodes( tre)[(length(tre$tip.label)+1), 1:(length(tre$tip.label)+1)], tre$tip.label)
# modify edge length to represent evolutionary distance with rate 1e-3
tre$edge.length <- tre$edge.length * 1e-3
# treedater: 
td <- dater( tre, sts =sts )
td$tre
td$meanRate
```

