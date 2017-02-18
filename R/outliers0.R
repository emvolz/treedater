
outlier.tips <- function(td, alpha = .01){
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
	qu.df <- data.frame( taxon = td$tip.label[ td$edge[tipEdges,2]]
	  , q = unname( p.adjust( p[tipEdges], method = 'fdr' ) )
	  , p  = unname( p[tipEdges] ) 
	  , loglik = unname( lls[tipEdges])
	  , rates = unname( td$omegas[tipEdges] )
	  , branch.length = td$edge.length[ tipEdges] 
	)
	qu.df <- qu.df[ order( qu.df$q), ]
	rownames(qu.df) <- qu.df$taxon
	qu.df <- qu.df[,-1] 
	qu.df1 <- qu.df[ qu.df$q <  alpha , ]
	
	print( qu.df1 )
	qu.df
}

