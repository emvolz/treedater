require(foreach)


parboot.treedater <- function( td , nreps = 100,  overrideTempConstraint=T, overrideClock=NULL, overrideSearchRoot=TRUE, overrideSeqLength = NULL, quiet=TRUE, normalApproxTMRCA=F, ncpu = 1, parallel_foreach = FALSE )
{
	if (quiet){
	cat( 'Running in quiet mode. To print progress, set quiet=FALSE.\n')
	}
	if (overrideSearchRoot){
		cat('NOTE: Running with overrideSearchRoot will speed up execution but may underestimate variance.\n')
	}
	if (overrideTempConstraint){
		cat('NOTE: Running with overrideTempConstraint will speed up execution but may underestimate variance.\n')
	}
	level <- .95
	alpha <- min(1, max(0, 1 - level ))

	.parboot.replicate <- function( k  = NA )
	{
		tre <- list( edge = td$edge, edge.length = td$edge.length, Nnode = td$Nnode, tip.label = td$tip.label)
			class(tre) <- 'phylo'
			blen <- pmax( 1e-12, td$edge.length )
			if( td$clock == 'strict' ) {
				# simulate poisson 
				tre$edge.length <- rpois(length(tre$edge.length), blen * td$meanRate * td$s ) / td$s
			} else {
				#ps <- pmin(1 - 1e-5, td$theta*blen  / ( 1+ td$theta * blen ) )
				ps <- pmin(1 - 1e-12, td$theta*blen  / ( 1+ td$theta * blen ) )
				tre$edge.length <- rnbinom( length(td$edge.length)
				 , td$r
				 , 1 - ps
				) / td$s
			}
			if (!td$intree_rooted & !overrideSearchRoot) tre <- unroot( tre )
			est <- NULL
			if (td$EST_SAMP_TIMES) est <- td$estimateSampleTimes
			tempConstraint <- td$temporalConstraints
			if ( overrideTempConstraint) tempConstraint <- FALSE
			clockstr <- td$clock
			if (!is.null( overrideClock)){
				if (is.na(overrideClock)) stop('overrideClock NA. Quitting.')
				if (!overrideClock %in% c('strict', 'relaxed') ){
					stop('overrideClock must be one of strict or relaxed')
				}
				clockstr <- overrideClock
			}
			strictClock <- ifelse( clockstr=='strict' , TRUE, FALSE )
			seqlen <- ifelse( is.null(overrideSeqLength) , td$s, overrideSeqLength )
			td2 <- tryCatch({dater(tre, td$sts, s= seqlen
			 , omega0 = NA
			 , minblen = td$minblen
			 , quiet = TRUE
			 , temporalConstraints = tempConstraint
			 , strictClock = strictClock
			 , estimateSampleTimes = est
			 , estimateSampleTimes_densities = td$estimateSampleTimes_densities
			 , numStartConditions = td$numStartConditions 
			 , meanRateLimits = td$meanRateLimits
			)}, error =function(e) NA )
			if (suppressWarnings( is.na( td2[1])) ) return (NA )
			if (!quiet){
				cat('\n #############################\n')
				cat( paste( '\n Replicate', k, 'complete \n' ))
				print( td2 )
			}
			td2
	}

	if (ncpu > 1)
	{
		if (parallel_foreach){
			tds <- foreach( k = 1:nreps, .packages=c('treedater') ) %dopar% .parboot.replicate(k)
		} else{
			tds <- parallel::mclapply( 1:nreps, function(k) .parboot.replicate(k) 
			, mc.cores = ncpu ) 
		}
	} else{
		tds <- foreach( k = 1:nreps, .packages=c('treedater') ) %do% .parboot.replicate( k)#td, overrideSearchRoot, overrideClock , quiet )
	}
	tds <- tds[!suppressWarnings ( sapply(tds, function(td) is.na(td[1]) ) )  ] 
	if (length(tds)==0) stop('All bootstrap replicate failed with error.')
	# output rate CIs, parameter CIs, trees
	meanRates <- sapply( tds, function(td) td$meanRate )
	cvs <- sapply( tds, function(td) td$coef_of_variation )
	tmrcas <- sapply( tds, function(td) td$timeOfMRCA )
	ttmrcas <- sapply( tds, function(td) td$timeTo )
	log_mr_sd <- sd(log(meanRates))
	log_tmrca_sd <- sd( log(ttmrcas ) )
	if (normalApproxTMRCA){
		timeOfMRCA_CI <- c( td$timeOfMRCA * exp(-log_tmrca_sd*1.96), td$timeOfMRCA * exp(log_tmrca_sd*1.96 ))
	} else {
		timeOfMRCA_CI <- quantile( tmrcas, p = c(.025, .975 ))
	}
	if (td$timeOf < timeOfMRCA_CI[1] | td$timeOf > timeOfMRCA_CI[2] ){
		warning( 'Parametric bootstrap CI does not cover the point estimate. Try non-parametric bootstrap *boot.treedater*. ')
	}
	rv <- list( 
		trees = tds
		, meanRates = meanRates
		, meanRate_CI = c( exp( log(td$meanRate) - log_mr_sd*1.96), exp(log(td$meanRate) + log_mr_sd*1.96 ))
		, coef_of_variation_CI = quantile( cvs, probs = c(alpha/2, 1 - alpha/2))
		, timeOfMRCA_CI = timeOfMRCA_CI
		, td = td
		, alpha = alpha
		, level = level
		, tmrcas = tmrcas
		, ttmrcas = ttmrcas
	)
	class(rv) <- 'boot.treedater'
	rv
}

print.parboot.treedater = print.boot.treedater <- function( x, ... )
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

plot.parboot.ltt = plot.boot.ltt <- function(pbtd, t0 = NA, res = 100, ... )
{
	t1 <- max( pbtd$td$sts, na.rm=T )
	if (is.na(t0)) t0 <- min( sapply( pbtd$trees, function(tr) tr$timeOf ) )
	times <- seq( t0, t1, l = res )
	cbind( times = times , t( sapply( times, function(t){
		c( pml = sum(pbtd$td$sts > t ) - sum( pbtd$td$Ti>t ) 
			, setNames(quantile( sapply( pbtd$trees, function(tre ) sum( tre$sts > t) - sum( tre$Ti > t ) )
			 , probs = c( .025, .5, .975 ) 
			), c('lb', 'median', 'ub') )
		)
	}))) -> ltt

	pl.df <- as.data.frame( ltt )
	p <- ggplot( pl.df ) + geom_ribbon( aes(x = times, ymin = lb, ymax = ub ), fill='blue', col = 'blue', alpha = .1 )
	p <- p + geom_path( aes(x = times, y = pml ))
	(p <- p + ylab( 'Lineages through time') + xlab('Time')  )
}



relaxed.clock.test <- function( ..., nreps=100, overrideTempConstraint=T )
{
	argnames <- names(list(...))
	if ( 'strictClock'  %in% argnames) stop('Can not prespecify clock type *strictClock* for relaxed.clock.test. Quitting.')
	td <- dater(..., strict=TRUE)
	pbtd <-  parboot.treedater( td , nreps = nreps,  overrideTempConstraint=overrideTempConstraint
	 , overrideClock = 'relaxed' )
	cvci_null <- pbtd$coef_of_variation_CI
	
	tdrc <- dater(..., strict=FALSE)
	clock <- 'strict'
	if ( tdrc$coef_of_variation > cvci_null[2] ){
		clock <- 'relaxed'
	}
	cat( paste( 'Best clock model: ', clock, '\n'))
	cat( paste( 'Null distribution of rate coefficient of variation:', paste(cvci_null, collapse=' '), '\n'))
	cat('Returning best treedater fit\n')
	list( strict_treedater = td
	 , relaxed_treedater = tdrc 
	 , clock = clock
	 , parboot_strict = pbtd
	 , nullHypothesis_coef_of_variation_CI = cvci_null
	)
}


##
boot.treedater <- function( td, tres, overrideTempConstraint=TRUE, searchRoot=1 , overrideClock=NULL, nreps = 100, quiet=TRUE, normalApproxTMRCA=F, ncpu = 1, parallel_foreach = FALSE )
{
	if (quiet){
		cat( 'Running in quiet mode. To print progress, set quiet=FALSE.\n')
	}
	if (overrideTempConstraint){
		cat('NOTE: Running with overrideTempConstraint will speed up execution but may underestimate variance.\n')
	}
	if (length(tres) < nreps){
		warning('Number of non-parametric bootstrap trees is less than number of bootstrap replicates requested. Results will under-estimate variance.')
	}
	level <- .95
	alpha <- min(1, max(0, 1 - level ))

	.boot.replicate <- function(  )
	{
		k <- sample.int( length(tres), size=1)
		tre <- tres[[k]]
		
		est <- NULL
		if (td$EST_SAMP_TIMES) est <- td$estimateSampleTimes
		tempConstraint <- td$temporalConstraints
		if ( overrideTempConstraint) tempConstraint <- FALSE
		clockstr <- td$clock
		if (!is.null( overrideClock)){
			if (is.na(overrideClock)) stop('overrideClock NA. Quitting.')
			if (!overrideClock %in% c('strict', 'relaxed') ){
				stop('overrideClock must be one of strict or relaxed')
			}
			clockstr <- overrideClock
		}
		strictClock <- ifelse( clockstr=='strict' , TRUE, FALSE )
		td2 <- tryCatch({ dater(tre, td$sts, s= td$s
		 , omega0 = NA
		 , minblen = td$minblen
		 , quiet = TRUE
		 , searchRoot = searchRoot
		 , temporalConstraints = tempConstraint
		 , strictClock = strictClock
		 , estimateSampleTimes = est
		 , estimateSampleTimes_densities = td$estimateSampleTimes_densities
		 , numStartConditions = td$numStartConditions 
		 , meanRateLimits = td$meanRateLimits
		)}, error =function(e) NA )
		if (suppressWarnings( is.na( td2[1])) ) return (NA )
			
		if (!quiet){
			cat('\n #############################\n')
			cat( paste( '\n Replicate', k, 'complete \n' ))
			print( td2 )
		}
		td2
	}

	if (ncpu > 1)
	{
		if (parallel_foreach){
			tds <- foreach( k = 1:nreps, .packages=c('treedater') ) %dopar% .boot.replicate()
		} else{
			tds <- parallel::mclapply( 1:nreps, function(k) .boot.replicate() 
			, mc.cores = ncpu ) 
		}
	} else{
		tds <- foreach( k = 1:nreps, .packages=c('treedater') ) %do% .boot.replicate( )
	}
	tds <- tds[!suppressWarnings ( sapply(tds, function(td) is.na(td[1]) ) )  ] 
	if (length(tds)==0) stop('All bootstrap replicate failed with error.')
	# output rate CIs, parameter CIs, trees
	meanRates <- sapply( tds, function(td) td$meanRate )
	cvs <- sapply( tds, function(td) td$coef_of_variation )
	tmrcas <- sapply( tds, function(td) td$timeOfMRCA )
	ttmrcas <- sapply( tds, function(td) td$timeTo )
	log_mr_sd <- sd(log(meanRates))
	log_tmrca_sd <- sd( log(ttmrcas ) )
	if (normalApproxTMRCA){
		timeOfMRCA_CI <- c( td$timeOfMRCA * exp(-log_tmrca_sd*1.96), td$timeOfMRCA * exp(log_tmrca_sd*1.96 ))
	} else {
		timeOfMRCA_CI <- quantile( tmrcas, p = c(.025, .975 ))
	}
	rv <- list( 
		trees = tds
		, meanRates = meanRates
		, meanRate_CI = c( exp( log(td$meanRate) - log_mr_sd*1.96), exp(log(td$meanRate) + log_mr_sd*1.96 ))
		, coef_of_variation_CI = quantile( cvs, probs = c(alpha/2, 1 - alpha/2))
		, timeOfMRCA_CI = timeOfMRCA_CI
		, td = td
		, alpha = alpha
		, level = level
		, tmrcas = tmrcas
		, ttmrcas = ttmrcas
	)
	class(rv) <- 'boot.treedater'
	rv
}
