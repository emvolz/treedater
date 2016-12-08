print.outlier.tips0 <- function(td, alpha = .1)
{
	w <- td$omegas
	rw <- range(w)
	if (isTRUE(all.equal(rw[1], rw[2]))){
		cat('Fit object does not show rate variation. Exiting.\n')
		return(NULL)
	}
	n <- length( td$tre$tip.label )
	tipEdges <- which( td$tre$edge[,2] <= n )
	qu <- sapply( 1:length(tipEdges), function(k){
		ww <- w[k] 
		blen <- max(td$minblen, td$tre$edge.length[k]  )
		min( pgamma( ww*td$s*blen, shape= td$r, scale = td$theta*blen ), 
		  1- pgamma( ww*td$s*blen, shape = td$r, scale = td$theta * blen )
		)
	})
	qu.df <- data.frame( taxon = td$tre$tip.label[ td$tre$edge[tipEdges,2]]
	  , p = unname( qu ) 
	)
	qu.df <- qu.df[ order( qu.df$p), ]
	qu.df1 <- qu.df[ qu.df$p < alpha, ]
	
	if (nrow(qu.df1)==0) cat( 'No outliers detected.') 
	print( qu.df1 )
	qu.df1
}

print.outlier.tips1 <- function(td, alpha = .1)
{
	n <- length( td$tre$tip.label )
	tipEdges <- which( td$tre$edge[,2] <= n )
	qu.df <- data.frame( taxon = td$tre$tip.label[ td$tre$edge[tipEdges,2]]
	  , loglik = unname( td$edge_lls[tipEdges]) 
	)
	qu.df <- qu.df[ order( qu.df$loglik), ]
	qu.df1 <- qu.df[ qu.df$loglik < quantile(qu.df$loglik, probs = alpha) , ]
	
	print( qu.df1 )
	qu.df
}

print.outlier.tips <- print.outlier.tips1
