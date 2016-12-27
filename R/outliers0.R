
outlier.tips <- function(td, alpha = .05){
	if (td$clock=='relaxed'){
		with(td,{
			blen <- pmax( minblen, tre$edge.length )
			ps <- pmin(1 - 1e-5, theta*blen / ( 1+ theta * blen ) )
			dnbinom( pmax(0, round(intree$edge.length*s))
			 , size= r, prob=1-ps,  log = T) 
		}) -> lls
		with(td,{
			blen <- pmax( minblen, tre$edge.length )
			ps <- pmin(1 - 1e-5, theta*blen / ( 1+ theta * blen ) )
			pnbinom( pmax(0, round(intree$edge.length*s))
			 , size= r, prob=1-ps) 
		}) -> p
	} else{
		with(td,{
			blen <- pmax( minblen, tre$edge.length )
			dpois( pmax(0, round(intree$edge.length*s))
			 , blen * meanRate *s,  log = T) 
		}) -> lls
		with(td,{
			blen <- pmax( minblen, tre$edge.length )
			ppois( pmax(0, round(intree$edge.length*s))
			 , blen * meanRate * s)
		}) -> p
	}
	p[ p > .5] <- 1 - p[ p > .5]
	n <- length( td$tre$tip.label )
	tipEdges <- which( td$tre$edge[,2] <= n )
	qu.df <- data.frame( taxon = td$tre$tip.label[ td$tre$edge[tipEdges,2]]
	  , loglik = unname( lls[tipEdges])
	  , p  = unname( p[tipEdges] ) 
	  , rates = unname( td$omegas[tipEdges] )
	  , branch.length = td$tre$edge.length[ tipEdges] 
	)
	qu.df <- qu.df[ order( qu.df$loglik), ]
	rownames(qu.df) <- qu.df$taxon
	qu.df <- qu.df[,-1] 
	qu.df1 <- qu.df[ qu.df$p <  alpha , ]
	
	print( qu.df1 )
	qu.df
}

