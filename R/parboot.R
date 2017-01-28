
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
	mr_sd <- sd(meanRates)
	tmrca_sd <- sd( tmrcas )
	rv <- list( 
		trees = tds
		, meanRates = meanRates
		, meanRate_CI = c( max(0, td$meanRate - mr_sd*1.96), td$meanRate + mr_sd*1.96 )
		#, meanRate_CI = quantile( meanRates, probs = c(alpha/2, 1-alpha/2 ))
		, coef_of_variation_CI = quantile( cvs, probs = c(alpha/2, 1 - alpha/2))
		, timeOfMRCA_CI = c( td$timeOfMRCA - tmrca_sd*1.96, td$timeOfMRCA + tmrca_sd*1.96 )
		#, timeOfMRCA_CI = quantile( tmrcas, probs = c(alpha/2, 1 - alpha/2))
		, td = td
		, alpha = alpha
		, level = level
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
