#' treedater fits a molecular clock to a phylogenetic tree and estimates evolutionary rates and times of common ancestry. 
#'
#' Additional functions are provided for detecting outlier lineages (possible sequencing or alignment error). A statistical test is available for choosing between strict and relaxed clock models. The calendar time of each sample must be specified (possibly with bounds of uncertainty) and the length of the sequences used to estimate the tree. treedater uses heuristic search to optimise the TMRCAs of a phylogeny and the substitution rate. An uncorrelated relaxed molecular clock accounts for rate variation between lineages of the phylogeny which is parameterised using a Gamma-Poisson mixture model.
#' @references
#' Volz, E. M., and S. D. W. Frost. "Scalable relaxed clock phylogenetic dating." Virus Evolution 3.2 (2017).
#' @import ape
#' @import limSolve
#' @import stats
#' @import utils
"_PACKAGE"
#> [1] "_PACKAGE"

