
outlier.lineages <- function(td, alpha = .01, type=c('tips','internal', 'all')){
	if (length(type)>1) type <- type[1]
	if ( !td$temporalConstraints ){
		stop('The outlier.tips function requires a treedater object fitted using temporalConstraints=TRUE. Quitting.')
	}
	if (td$clock=='relaxed'){
		with(td,{
			blen <- pmax( minblen, edge.length )
			ps <- pmin(1 - 1e-5, theta*blen / ( 1+ theta * blen ) )
			dnbinom( pmax(0, round(intree$edge.length*s))
			 , size= r, prob=1-ps,  log = T) 
		}) -> lls
		with(td,{
			blen <- pmax( minblen, edge.length )
			ps <- pmin(1 - 1e-5, theta*blen / ( 1+ theta * blen ) )
			pnbinom( pmax(0, round(intree$edge.length*s))
			 , size= r, prob=1-ps) 
		}) -> p
	} else{
		with(td,{
			blen <- pmax( minblen, edge.length )
			dpois( pmax(0, round(intree$edge.length*s))
			 , blen * meanRate *s,  log = T) 
		}) -> lls
		with(td,{
			blen <- pmax( minblen, edge.length )
			ppois( pmax(0, round(intree$edge.length*s))
			 , blen * meanRate * s)
		}) -> p
	}
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



outlier.tips <- function( td, alpha = .01){
	outlier.lineages(td, alpha = .01)
}
