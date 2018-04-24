## ------------------------------------------------------------------------
require(treedater)
(tre <- ape::read.tree( system.file( 'extdata', 'flu_h3n2_final_small.treefile', package='treedater') ))

## ------------------------------------------------------------------------
seqlen <- 1698 # the length of the HA sequences used to reconstruct the phylogeny

## ------------------------------------------------------------------------
sts <- sampleYearsFromLabels( tre$tip.label, delimiter='_' )
head(sts)

## ------------------------------------------------------------------------
hist( sts , main = 'Time of sequence sampling') 

