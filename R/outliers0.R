#~ Treedater: fast relaxed molecular clock dating 
#~     Copyright (C) 2018  Erik Volz
#~     This program is free software: you can redistribute it and/or modify
#~     it under the terms of the GNU General Public License as published by
#~     the Free Software Foundation, either version 3 of the License, or
#~     (at your option) any later version.


#' Detect lineages with unusually large evolutionary divergence under the fitted treedater model
#'
#' Outliers are detected using the *stats::p.adjust* function and the 'fdr' function. The test requires that *dater* was used with the temporalConstraints=TRUE.
#'
#' @param td A fitted treedater object
#' @param alpha The tail probability used for classifying lineages as outliers
#' @param type Should outliers be detected on tip lineages, interal lineages, or all lineages? 
#' @return A data frame summarizing for each lineage the p values, adjusted p values ('q'), likelihood, rates, and branch lengths. 
#' @seealso 
#' dater
#' outlier.tips
#' @export
outlierLineages <- function(td, alpha = .05, type=c('tips','internal', 'all')){
	if (length(type)>1) type <- type[1]
	if ( !td$temporalConstraints ){
		stop('The outlier.tips function requires a treedater object fitted using temporalConstraints=TRUE. Quitting.')
	}
	lls <- td$edge_lls
	p <- td$edge.p
	p[ p > .5] <- 1 - p[ p > .5]
	p <- p * 2 
	n <- length( td$tip.label )
	tipEdges <- which( td$edge[,2] <= n )
	intEdges <- setdiff( 1:nrow(td$edge), tipEdges )
	if ( type == 'tips'){
		qu.df <- data.frame( taxon = td$tip.label[ td$edge[tipEdges,2]]
		  , q = unname( p.adjust( p[tipEdges], method = 'fdr' ) )
		  , p  = unname( p[tipEdges] ) 
		  , loglik = unname( lls[tipEdges])
		  , rates = unname( td$omegas[tipEdges] )
		  , branch.length = td$edge.length[ tipEdges] 
		)
		rownames(qu.df) <- qu.df$taxon
	} else if ( type == 'internal'){
		qu.df <- data.frame( 
		   q = unname( p.adjust( p[intEdges], method = 'fdr' ) )
		  , p  = unname( p[intEdges] ) 
		  , loglik = unname( lls[intEdges])
		  , rates = unname( td$omegas[intEdges] )
		  , branch.length = td$edge.length[ intEdges ] 
		  , edgeIndex = intEdges
		)
		#rownames(qu.df) <- qu.df$intEdges
	} else{
		qu.df <- data.frame(
		   q = unname( p.adjust( p, method = 'fdr' ) )
		  , p  = unname( p ) 
		  , loglik = unname( lls)
		  , rates = unname( td$omegas )
		  , branch.length = td$edge.length 
		  , edgeIndex = 1:nrow(td$edge)
		)
	}
	qu.df <- qu.df[ order( qu.df$q), ]
	
	qu.df1 <- qu.df[ qu.df$q <  alpha , ]
	
	print( qu.df1 )
	qu.df
}

#' Detect terminal lineages with unusually large evolutionary divergence under the fitted treedater model
#'
#' This is a convient wrapper of the *outlier.lineages*
#'
#' @param td A fitted treedater object
#' @param alpha The tail probability used for classifying lineages as outliers
#' @return A data frame summarizing for each lineage the p values, adjusted p values ('q'), likelihood, rates, and branch lengths. 
#' @seealso 
#' dater
#' outlier.lineages
#' @export
outlierTips <- function( td, alpha = .05){
	outlierLineages(td, alpha = alpha)
}
