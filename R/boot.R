#~ Treedater: fast relaxed molecular clock dating 
#~     Copyright (C) 2018  Erik Volz
#~     This program is free software: you can redistribute it and/or modify
#~     it under the terms of the GNU General Public License as published by
#~     the Free Software Foundation, either version 3 of the License, or
#~     (at your option) any later version.


#' Estimate of confidence intervals using parametric bootstrap for molecular clock dating.
#'
#'     This function simulates phylogenies with branch lengths in units
#'      of substitutions per site. Simulations are based on a fitted
#'      treedater object which provides parameters of the molecular clock
#'      model. The treedater method is applied to each simulated tree
#'      providing a Monte Carlo estimate of variance in rates and dates.
#'
#'      If the original treedater fit estimated the root position, root
#'      position will also be estimated for each simulation, so the
#'      returned trees may have different root positions. Some replicates
#'      may converge to a strict clock or a relaxed clock, so the
#'      parameter estimates in each replicate may not be directly
#'      comparable. It is possible to compute confidence intervals for the
#'      times of particular nodes or for estimated sample times by
#'      inspecting the output from each fitted treedater object, which is
#'      contained in the $trees attribute.
#'
#' @param td A fitted treedater object 
#' @param nreps Integer number of simulations to be carried out 
#' @param ncpu Number of threads to use for parallel computation. Recommended.
#' @param overrideTempConstraint If TRUE (default) will not enforce positive branch lengths in simualtion replicates. Will speed up execution. 
#' @param overrideClock May be 'strict' or 'relaxed' in which case will force simulations to fit the corresponding model. If ommitted, will inherit the clock model from td
#' @param overrideSearchRoot If TRUE, will re-use root position from input treedater tree. Otherwise may re-estimate root position in simulations
#' @param overrideSeqLength Optional sequence length to use in simulations
#' @param quiet If TRUE will minimize output printed to screen
#' @param normalApproxTMRCA If TRUE will use estimate standard deviation from simulation replicates and report confidence interval based on normal distribution
#' @param parallel_foreach If TRUE will use the foreach package for parallelization. May work better on HPC systems. 
#'
#' @return 
#' A list with elements 
#' \itemize{
#' \item trees: The fitted treedater objects corresponding to each simulation
#' \item meanRates: Vector of estimated rates for each simulation
#' \item meanRate_CI: Confidence interval for substitution rate
#' \item coef_of_variation_CI: Confidence interval for rate variation
#' \item timeOfMRCA_CI: Confidence interval for time of common ancestor
#' }
#'
#' @seealso
#' dater
#' boot
#'
#' @author Erik M Volz <erik.volz@gmail.com>
#'
#' @examples
#' # make a random tree
#' tre <- ape::rtree(25)
#' # simulate sample times based on distance from root to tip:
#' sts <- setNames( ape::node.depth.edgelength( tre )[1:ape::Ntip(tre)], tre$tip.label)
#' # modify edge length to represent evolutionary distance with rate 1e-3:
#' tre$edge.length <- tre$edge.length * 1e-3
#' # treedater: 
#' td <- dater( tre, sts =sts )
#' # parametric bootstrap: 
#' pb <- parboot( td, nreps=25 )
#' # plot lineages through time
#' plot( pb )
#'
#' @export 
parboot <- function( td , nreps = 100, ncpu = 1,  overrideTempConstraint=TRUE, overrideClock=NULL, overrideSearchRoot=TRUE, overrideSeqLength = NULL, quiet=TRUE, normalApproxTMRCA=FALSE, parallel_foreach = FALSE )
{
	k = 0 # resolve NOTE about 'visible binding for global variable'
	if (quiet){
	cat( 'Running in quiet mode. To print progress, set quiet=FALSE.\n')
	}
	if (overrideSearchRoot){
		cat('NOTE: Running with overrideSearchRoot will speed up execution but may underestimate variance.\n')
	}
	if (overrideTempConstraint){
		cat('NOTE: Running with overrideTempConstraint will speed up execution but may underestimate variance. Bootstrap tree replicates may have negative branch lengths.\n')
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
				if (is.na(overrideClock)) stop('overrideClock NA. Stopping.')
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
			success = requireNamespace('foreach', quietly=TRUE)
			if (!success) stop('*parallel_foreach* requires the `foreach` package. Stopping.')
			`%dopar%` <- foreach::`%dopar%`
			tds <- foreach::foreach( k = 1:nreps, .packages=c('treedater') ) %dopar% {
				capture.output( { pbrk <- .parboot.replicate(k) })
				pbrk 
			}
		} else{
			tds <- parallel::mclapply( 1:nreps, function(k) {
				capture.output( { pbrk <- .parboot.replicate(k) })
				pbrk
			}
			, mc.cores = ncpu ) 
		}
	} else{
		tds <- lapply( 1:nreps, function(k) {
			capture.output( { pbrk <- .parboot.replicate(k) })
			pbrk
		})
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
	class(rv) <- 'bootTreedater'
	rv
}


#' Estimate of confidence intervals of molecular clock parameters with user-supplied set of bootstrap trees
#'
#'      If the original treedater fit estimated the root position, root
#'      position will also be estimated for each simulation, so the
#'      returned trees may have different root positions. Some replicates
#'      may converge to a strict clock or a relaxed clock, so the
#'      parameter estimates in each replicate may not be directly
#'      comparable. It is possible to compute confidence intervals for the
#'      times of particular nodes or for estimated sample times by
#'      inspecting the output from each fitted treedater object, which is
#'      contained in the $trees attribute.
#'
#' @param td A fitted treedater object 
#' @param tres A list or multiPhylo with bootstrap trees with branches in units of substitutions per site 
#' @param ncpu Number of threads to use for parallel computation. Recommended.
#' @param searchRoot See *dater*
#' @param overrideTempConstraint If TRUE (default) will not enforce positive branch lengths in simualtion replicates. Will speed up execution. 
#' @param overrideClock May be 'strict' or 'relaxed' in which case will force simulations to fit the corresponding model. If ommitted, will inherit the clock model from td
#' @param quiet If TRUE will minimize output printed to screen
#' @param normalApproxTMRCA If TRUE will use estimate standard deviation from simulation replicates and report confidence interval based on normal distribution
#' @param parallel_foreach If TRUE will use the foreach package for parallelization. May work better on HPC systems. 
#'
#' @return 
#' A list with elements 
#' \itemize{
#' \item trees: The fitted treedater objects corresponding to each simulation
#' \item meanRates: Vector of estimated rates for each simulation
#' \item meanRate_CI: Confidence interval for substitution rate
#' \item coef_of_variation_CI: Confidence interval for rate variation
#' \item timeOfMRCA_CI: Confidence interval for time of common ancestor
#' }
#'
#' @seealso
#' dater
#' parboot
#'
#' @author Erik M Volz <erik.volz@gmail.com>
#'
#' @examples 
#' # simulate a tree 
#' tre <- ape::rtree(25)
#' # sample times based on distance from root to tip:
#' sts <- setNames( ape::node.depth.edgelength( tre )[1:ape::Ntip(tre)], tre$tip.label)
#' # make a list of trees that simulate outcome of bootstrap using nonparametric phylogeny estimation
#' # also modify edge length to represent evolutionary distance with rate 1e-3:
#' bootTrees <- lapply( 1:25, function(i) {
#' 	.tre <- tre
#' 	.tre$edge.length <- tre$edge.length * pmax(rnorm( length(tre$edge.length), 1e-3, 1e-4 ), 0 )
#' 	.tre 
#' })
#' tre$edge.length <- tre$edge.length * 1e-3
#' # run treedater
#' td <- dater( tre, sts  )
#' # bootstrap: 
#' ( tdboot <- boot( td, bootTrees ) )
#' # plot lineages through time :
#' plot( tdboot )
#'
#'
#' @export 
boot <- function( td, tres,  ncpu = 1, searchRoot=1 , overrideTempConstraint=TRUE,  overrideClock=NULL, quiet=TRUE, normalApproxTMRCA=FALSE, parallel_foreach = FALSE )
{
	k = 0 # resolve NOTE about 'visible binding for global variable'
	nreps <- length(tres )
	if (quiet){
		cat( 'Running in quiet mode. To print progress, set quiet=FALSE.\n')
	}
	if (overrideTempConstraint){
		cat('NOTE: Running with overrideTempConstraint will speed up execution but may underestimate variance. Bootstrap tree replicates may have negative branch lengths.\n')
	}
	if (length(tres) < nreps){
		warning('Number of non-parametric bootstrap trees is less than number of bootstrap replicates requested. Results will under-estimate variance.')
	}
	level <- .95
	alpha <- min(1, max(0, 1 - level ))

	.boot.replicate <- function( k = NULL )
	{
		if (is.null(k)) k <- sample.int( length(tres), size=1)
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
			success = requireNamespace('foreach', quietly=TRUE)
			if (!success) stop('*parallel_foreach* requires the `foreach` package. Stopping.')
			`%dopar%` <- foreach::`%dopar%`
			tds <- foreach::foreach( k = 1:nreps, .packages=c('treedater') ) %dopar% {
				 capture.output( { brk <- .boot.replicate(k) })
				 brk 
			}
		} else{
			tds <- parallel::mclapply( 1:nreps, function(k) {
				capture.output( { brk <- .boot.replicate(k) })
				brk 
			}
			, mc.cores = ncpu ) 
		}
	} else{
		tds <- lapply( 1:nreps, function(k) {
			capture.output( { brk <- .boot.replicate(k) })
			brk 
		})
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
	class(rv) <- 'bootTreedater'
	rv
}


#'Use parametric bootstrap to test if relaxed clock offers improved fit to data.
#'
#'      This function simulates phylogenies with branch lengths in units
#'      of substitutions per site. Simulations are based on a fitted
#'      treedater object which provides parameters of the molecular clock
#'      model. The coefficient of variation of rates is estimated using a
#'      relaxed clock model applied to strict clock simulations. Estimates
#'      of the CV is then compared to the null distribution provided by
#'      simulations.
#'
#'      This function will print the optimal clock model
#'      and the distribution of the coefficient of variation statistic under the null hypothesis (strict
#'      clock). Parameters passed to this function should be the same as when calling *dater*.
#' 
#' @param ... arguments passed to *dater*
#' @param nreps Integer number of simulations
#' @param overrideTempConstraint see *parboot*
#' @param ncpu Number of threads to use for parallel computation. Recommended.
#' 
#' @return  A list with elements:
#' \itemize{
#' \item strict_treedater: A dater object under a strict clock
#' \item relaxed_treedater: A dater object under a relaxed clock
#' \item clock: The favoured clock model 
#' \item parboot: Result of call to *parboot* using fitted treedater and forcing a relaxed clock
#' \item nullHypothesis_coef_of_variation_CI: The null hypothesis CV
#' }
#'
#' @author Erik M Volz <erik.volz@gmail.com>
#'
#' @examples 
#' # simulate a tree 
#' tre <- ape::rtree(25)
#' # sample times based on distance from root to tip:
#' sts <- setNames( ape::node.depth.edgelength( tre )[1:ape::Ntip(tre)], tre$tip.label)
#' # modify edge length to represent evolutionary distance with rate 1e-3:
#' tre$edge.length <- tre$edge.length * 1e-3
#' relaxedClockTest( tre, sts, nreps=25)
#'
#'
#' @export 
relaxedClockTest <- function( ..., nreps=100, overrideTempConstraint=T , ncpu =1 )
{
	argnames <- names(list(...))
	if ( 'strictClock'  %in% argnames) stop('Can not prespecify clock type *strictClock* for relaxed.clock.test. Quitting.')
	td <- dater(..., strictClock=TRUE)
	pbtd <-  parboot( td , nreps = nreps,  overrideTempConstraint=overrideTempConstraint
	 , overrideClock = 'relaxed' , ncpu = ncpu )
	cvci_null <- pbtd$coef_of_variation_CI
	
	tdrc <- dater(..., strictClock=FALSE)
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
	 , parboot = pbtd
	 , nullHypothesis_coef_of_variation_CI = cvci_null
	)
}

##

#' @export 
print.bootTreedater = print.boot.treedater <- function( x, ... )
{
	stopifnot(inherits(x, "bootTreedater"))
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

#' Plots lineages through time and confidence intervals estimated by bootstrap. 
#'
#' @param x A bootTreedater object produced by *parboot* or *boot*
#' @param t0 The lower bound of the time axis to show
#' @param res The number of time points on the time axis
#' @param ggplot If TRUE, will return a plot made with the ggplot2 package
#' @param cumulative If TRUE, will show only decreasing lineages through time
#' @param ... Additional arg's are passed to *ggplot* or *plot*
#' @export 
plot.bootTreedater <- function(x, t0 = NA, res = 100, ggplot=FALSE, cumulative=FALSE, ... )
{
	pbtd = x
	stopifnot(inherits(pbtd, "bootTreedater"))
	t1 <- max( pbtd$td$sts, na.rm=T )
	if (is.na(t0)) t0 <- min( sapply( pbtd$trees, function(tr) tr$timeOf ) )
	times <- seq( t0, t1, l = res )
	if (cumulative)
	{
		cbind( times = times , t( sapply( times, function(t){
			c( pml = length(pbtd$td$sts) - sum( pbtd$td$Ti>t ) 
				, setNames(quantile( sapply( pbtd$trees, function(tre ) length(pbtd$td$sts) - sum( tre$Ti > t ) )
				 , probs = c( .025, .5, .975 ) 
				), c('lb', 'median', 'ub') )
			)
		}))) -> ltt
	} else{
		cbind( times = times , t( sapply( times, function(t){
			c( pml = sum(pbtd$td$sts > t ) - sum( pbtd$td$Ti>t ) 
				, setNames(quantile( sapply( pbtd$trees, function(tre ) sum( tre$sts > t) - sum( tre$Ti > t ) )
				 , probs = c( .025, .5, .975 ) 
				), c('lb', 'median', 'ub') )
			)
		}))) -> ltt
	}
	# resolve NOTE about 'no visible binding for global variables'
	pml = NULL
	ub = NULL 
	lb = NULL
	
	pl.df <- as.data.frame( ltt )
	if (ggplot){
		success = requireNamespace('ggplot2', quietly=TRUE)
		if (!success) stop('*ggplot* option requires that the `ggplot2` package is installed. Stopping.')
		p <- ggplot2::ggplot( pl.df, ... ) + ggplot2::geom_ribbon( ggplot2::aes(x = times, ymin = lb, ymax = ub ), fill='blue', col = 'blue', alpha = .1 )
		p <- p + ggplot2::geom_path( ggplot2::aes(x = times, y = pml ))
		return (p <- p + ggplot2::ylab( 'Lineages through time') + ggplot2::xlab('Time')  )
	}
	with( pl.df ,{
		graphics::plot( times, lb, type = 'l', lty = 3, lwd = 1, xlab = 'Time', ylab= 'Lineages through time'#, main='' 
		  , ylim = c(0, max(ub)+1)  , ...)
		graphics::lines( times, ub, lty = 3, lwd = 1)
		graphics::lines( times, pml, lwd = 2) 
	})
}

