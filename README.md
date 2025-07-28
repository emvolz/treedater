# treedater
`treedater` fits a strict or relaxed molecular clock to a phylogenetic tree and estimates evolutionary rates and times of common ancestry. The calendar time of each sample must be specified (possibly with bounds of uncertainty) and the length of the sequences used to estimate the tree. 

`treedater` uses heuristic search to optimise the TMRCAs of a phylogeny and the substitution rate. 
An uncorrelated relaxed molecular clock accounts for rate variation between lineages of the phylogeny which is parameterised using a Gamma-Poisson mixture model.

To cite:
* E.M. Volz and Frost, S.D.W. (2017) [Scalable relaxed clock phylogenetic dating](https://doi.org/10.1093/ve/vex025). Virus Evolution.

## Installation
You can install the latest development version from github using the `devtools` package:
```
library(devtools)
install_github( 'emvolz/treedater')
```

## Basic usage

```
dater( tre, sts, s, omega0)
```

where 
* `tre` is an `ape::phylo` phylogeny, 
* `sts` is a named vector of sample times for each tip in `tre`
* `s` is the length of the genetic sequences used to estimate `tre`
* `omega0` is an initial guess of the substitution rate (can be omitted)

For a detailed introduction to features available in `treedater`, see the vignette on analysis of Influenza H3N2: `vignette('h3n2')`. 

## Command line

You can also use treedater from the command line without starting R using the `tdcl` script: 
```
./tdcl -h
Usage: ./tdcl [-[-help|h] [<logical>]] [-[-treefn|t] <character>] [-[-samplefn|s] <character>] [-[-sequenceLength|l] <double>] [-[-output|o] [<character>]]

-t <file> : file name of tree in newick format  
-s  <file> : should be a comma-separated-value file with sample times in format <taxon-id,sample-time> and no header
-l <length> :  the integer length of sequences in alignment used to construct the tree 
-o <file>: name of file for saving output 
```
Note that you may need to modify the first line of the `tdcl` script with the correct path to `Rscript` or `littler`.



