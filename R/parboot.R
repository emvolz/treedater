
parboot.treedater <- function( td , nreps = 100,  overrideTempConstraint=T )
{
	level <- .95
	alpha <- min(1, max(0, 1 - level ))
	lapply( 1:nreps, function(k) 
	{
		tre <- list( edge = td$edge, edge.length = td$edge.length, Nnode = td$Nnode, tip.label = td$tip.label)
		class(tre) <- 'phylo'
		blen <- pmax( td$minblen, td$edge.length )
		if( td$clock == 'strict' ) {
			# simulate poisson 
			tre$edge.length <- rpois(length(tre$edge.length), blen * td$meanRate * td$s ) / td$s
		} else {
			ps <- pmin(1 - 1e-5, td$theta*blen  / ( 1+ td$theta * blen ) )
			tre$edge.length <- rnbinom( length(tre$edge.length)
			 , td$r
			 , 1 - ps
			) / td$s
		}
		if (!td$intree_rooted) tre <- unroot( tre )
		est <- NULL
		if (td$EST_SAMP_TIMES) est <- td$estimateSampleTimes
		tempConstraint <- td$temporalConstraints
		if ( overrideTempConstraint) tempConstraint <- FALSE
		td <- tryCatch({ dater(tre, td$sts
		 , omega0 = NA#td$meanRate
		 , minblen = td$minblen
		 , quiet = TRUE
		 , temporalConstraints = tempConstraint
		 , strictClock = ifelse( td$clock=='strict' , TRUE, FALSE )
		 , estimateSampleTimes = est
		 , estimateSampleTimes_densities = td$estimateSampleTimes_densities
		 , numStartConditions = td$numStartConditions 
		)}, error = function(e) NA)
		if (suppressWarnings( is.na( td[1])) ) return (NA )
		cat('\n #############################\n')
		cat( paste( '\n Replicate', k, 'complete \n' ))
		print( td )
		td
	}) -> tds
	tds <- tds[!suppressWarnings ( sapply(tds, function(td) is.na(td[1]) ) )  ] 
	# output rate CIs, parameter CIs, trees
	meanRates <- sapply( tds, function(td) td$meanRate )
	cvs <- sapply( tds, function(td) td$coef_of_variation )
	tmrcas <- sapply( tds, function(td) td$timeOfMRCA )
	ttmrcas <- sapply( tds, function(td) td$timeTo )
	log_mr_sd <- sd(log(meanRates))
	log_tmrca_sd <- sd( log(ttmrcas ) )
	rv <- list( 
		trees = tds
		, meanRates = meanRates
		, meanRate_CI = c( exp( log(td$meanRate) - log_mr_sd*1.96), exp(log(td$meanRate) + log_mr_sd*1.96 ))
		, coef_of_variation_CI = quantile( cvs, probs = c(alpha/2, 1 - alpha/2))
		, timeOfMRCA_CI = c( td$timeOfMRCA * exp(-log_tmrca_sd*1.96), td$timeOfMRCA * exp(log_tmrca_sd*1.96 ))
		, td = td
		, alpha = alpha
		, level = level
		, tmrcas = tmrcas
		, ttmrcas = ttmrcas
	)
	class(rv) <- 'parboot.treedater'
	rv
}

print.parboot.treedater <- function( x, ... )
{
rns <- c( 'Time of common ancestor' 
  , 'Mean substitution rate'
)
cns <- c( 'pseudo ML'
 , paste( round(100*x$alpha/2, digits=3), '%')
 , paste( round(100*(1-x$alpha/2), digits=3), '%')
)
o <- data.frame( pseudoML= c(x$td$timeOfMRCA, x$td$meanRate)
 , c( x$timeOfMRCA_CI[1], x$meanRate_CI[1])  
 , c( x$timeOfMRCA_CI[2], x$meanRate_CI[2] )  
)
	rownames(o) <- rns
	colnames(o) <- cns
	print(o)
	cat( '\n For more detailed output, $trees provides a list of each fit to each simulation \n')
	invisible(x)
}

plot.parboot.ltt <- function(pbtd, t0 = NA, res = 100, ... )
{
	t1 <- max( pbtd$sts )
	if (is.na(t0)) t0 <- min( sapply( pbtd$trees), function(tr) tr$timeOf )
	times <- seq( t0, t1, l = res )
	cbind( times = times , t( sapply( times, function(t){
		c( pml = sum(pbtd$td$sts > t ) - sum( pbtd$td$Ti>t ) 
			, setNames(quantile( sapply( pbtd$trees, function(tre ) sum( tre$sts > t) - sum( tre$Ti > t ) )
			 , probs = c( .025, .5, .975 ) 
			), c('lb', 'median', 'ub') )
		)
	}))) -> ltt

	pl.df <- as.data.frame( ltt )
	p <- ggplot( pl.df, ... ) + geom_ribbon( aes(x = times, ymin = lb, ymax = ub ), fill='blue', col = 'blue', alpha = .1 )
	p <- p + geom_path( aes(x = times, y = pml ))
	(p <- p + ylab( 'Lineages through time') + xlab('Time' )  )
}
