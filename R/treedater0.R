#~ Treedater: fast relaxed molecular clock dating 
#~     Copyright (C) 2018  Erik Volz
#~     This program is free software: you can redistribute it and/or modify
#~     it under the terms of the GNU General Public License as published by
#~     the Free Software Foundation, either version 3 of the License, or
#~     (at your option) any later version.


#' Compute a vector of numeric sample times from labels in a sequence aligment or phylogeny
#' 
#' @param tips A character vector supplying the name of each sample 
#' @param dateFormat The format of the sample date. See ?Date for more information
#' @param delimiter Character(s) which separate data in each label
#' @param index Integer position of the date string in each label with respect to *delimiter*
#' @param regex A regular expression for finding the date substring. Should not be used with *delimiter* or *index*
#' @return Numeric vector with sample time in decimal format. 
#' @examples 
#' ## A couple of labels for Ebola virus sequences: 
#' sampleYearsFromLabels( c('EBOV|AA000000|EM104|SierraLeone_EM|2014-06-02'
#'                        , 'EBOV|AA000000|G3713|SierraLeone_G|2014-06-09')
#'	, delimiter='|' )
#' ## Equivalently: 
#' sampleYearsFromLabels( c('EBOV|AA000000|EM104|SierraLeone_EM|2014-06-02'
#'                        , 'EBOV|AA000000|G3713|SierraLeone_G|2014-06-09')
#'  , regex='[0-9]+-[0-9]+-[0-9]+')
#'
#' @export 
sampleYearsFromLabels <- function(tips, dateFormat='%Y-%m-%d'
 , delimiter=NULL #'|'
 , index=NULL
 , regex=NULL
)
{
	success = requireNamespace('lubridate', quietly=TRUE)
	if (!success) stop(' The `sampleYearsFromLabels` function requires the `lubridate` package. Please install `lubridate`. Stopping.')
	
	#units <- match.arg(units)
	if (!is.null(delimiter) & !is.null(regex))
	 stop('At most one of the options *delimiter* or *regex* should be supplied for determining times of sampled lineages. The *index* parameter may be supplied with delimiter, or if omitted, will assume that the last field in the tip label corresponds to date.')
	
	if(!is.null(delimiter)){
	  if(is.null(index)){
		tipdates <- sapply( strsplit( tips, delimiter, fixed=TRUE) , function(x) tail(x,1) )
	  } else {
	    tipdates <- sapply( strsplit( tips, delimiter, fixed=TRUE) , function(x) x[index] )
	  }
	} else if( !is.null(regex)){
		tipdates <- regmatches( tips, regexpr( regex, tips ))
		if (length( tipdates)!= length(tips)) 
		  stop('Error with regex: number of matches does not equal number of tips.')
	}
	
	tipdates <- as.Date( tipdates, format=dateFormat )
	sts <- lubridate::decimal_date( tipdates )
	setNames( sts, tips )
}

.make.tree.data <- function( tre, sampleTimes, s, cc)
{
	n <- length( tre$tip.label)
	
	tipEdges <- which( tre$edge[,2] <= n)
	i_tip_edge2label <- tre$tip.label[ tre$edge[tipEdges,2] ]
	sts <- sampleTimes[i_tip_edge2label]
	
	# daughters, parent
	daughters <- matrix( NA, nrow = n + n-1, ncol = 2)
	parent <- rep(NA, n + n - 1)
	for (k in 1:nrow(daughters)){
		x <- tre$edge[which(tre$edge[,1] == k),2]
		if (length(x) > 0){
			daughters[k, ] <- x
			for (u in x){
				if (!is.na(u)) parent[u] <- k
			}
		}
	}
	
	B <- tre$edge.length
	A <- matrix(0, nrow = length(B), ncol = n-1)
	
	#B[tipEdges] <- B[tipEdges] -  unname( omega * sts )
	A[ cbind(1:nrow(tre$edge), tre$edge[,1]-n) ] <- -1
	internalEdges <- setdiff( 1:length(B), tipEdges )
	A[ cbind(internalEdges, tre$edge[internalEdges,2] - n) ] <- 1
	
	# constraints(optional)
	Ain <-  matrix(0, nrow = length(B), ncol = n-1)
	Ain[ cbind(1:nrow(tre$edge), tre$edge[,1]-n) ] <- 1 # parent to -1
	Ain[ cbind(internalEdges, tre$edge[internalEdges,2] - n) ] <- -1 # dgtr to +1
	bin <- rep(0, length(B))
	bin[tipEdges] <- sts # terminal edges to -sample time #...
	
	W <- abs( 1/( (tre$edge.length + cc / s)/s) ) 
	#W <- abs( 1/( pmax(.001,tre$edge.length) )/s)
	
	list( A0 = A, B0 = B, W0 = W, n = n, tipEdges=tipEdges
	 , i_tip_edge2label = i_tip_edge2label
	 , sts2 = sts  # in order of tipEdges
	 , sts1 = sampleTimes[ tre$tip.label] # in order of tip.label
	 , sts = sampleTimes
	 , s = s, cc = cc, tre = tre
	 , daughters = daughters
	 , parent = parent
	 , Ain = Ain
	 , bin = bin
	)
}

.Ti2blen <- function(Ti, td ){
	Ti <- c( td$sts[td$tre$tip.label], Ti)
	elmat <- cbind( Ti[td$tre$edge[,1]], Ti[td$tre$edge[,2]])
	pmax(td$minblen,-elmat[,1]  + elmat[,2] )
}

.optim.r.gammatheta.nbinom0 <- function(  Ti, r0, gammatheta0, td, lnd.mean.rate.prior)
{	
	blen <- .Ti2blen( Ti, td )
	
	#NOTE relative to wikipedia page on NB:
	#size = r
	#1-prob = p
	of <- function(x)
	{
		r <- max(1e-3, exp( x['lnr'] ) )
		gammatheta <- exp( x['lngammatheta'] )
		if (is.infinite(r) | is.infinite(gammatheta)) return(Inf)
		ps <- pmin(1 - 1e-5, gammatheta*blen / ( 1+ gammatheta * blen ) )
		ov <- -sum( dnbinom( pmax(0, round(td$tre$edge.length*td$s))
		  , size= r, prob=1-ps,  log = T) )
		mr <-  r * gammatheta / td$s # Note this value will differ slightly from output of .mean.rate
		ov <- ov - unname(lnd.mean.rate.prior( mr ))
		ov 
	}
	x0 <- c( lnr = unname(log(r0)), lngammatheta = unname( log(gammatheta0)))
	#o <- optim( par = x0, fn = of, method = 'BFGS' )
	if ( is.infinite( of( x0 )) ) stop( 'Can not optimize rate parameters from initial conditions. Try adjusting *meanRateLimits*.' )
	o <- optim( par = x0, fn = of)
	r <- unname( exp( o$par['lnr'] ))
	gammatheta <- unname( exp( o$par['lngammatheta'] ))
	list( r = r, gammatheta=gammatheta, ll = -o$value)
}

.optim.omega.poisson0 <- function(Ti, omega0, td, lnd.mean.rate.prior, meanRateLimits)
{	
	blen <- .Ti2blen( Ti, td )
	
	of <- function(omega)
	{
		-sum( dpois( pmax(0, round(td$tre$edge.length*td$s)), td$s * blen * omega ,  log = T) ) - unname(lnd.mean.rate.prior( omega ))
	}
	o <- optimise(  of, lower = max( omega0 / 10, meanRateLimits[1] )
	 , upper = min( omega0 * 10, meanRateLimits[2] ))
	list( omega = unname( o$minimum), ll = -unname(o$objective) )
}



.optim.Ti0 <- function( omegas, td , scale_var_by_rate = FALSE){
		A <- omegas * td$A0 
		B <- td$B0
		B[td$tipEdges] <- td$B0[td$tipEdges] -  unname( omegas[td$tipEdges] * td$sts2 )
		#solve( t(A) %*% A ) %*% t(A) %*% B
		if (scale_var_by_rate){
			rv <- ( coef( lm ( B ~ A -1 , weights = td$W/omegas) ) )
		} else{
			rv <- ( coef( lm ( B ~ A -1 , weights = td$W) ) ) 
		}
	if (any(is.na(rv))){
		warning('Numerical error when performing least squares optimisation. Values are approximate. Try adjusting minimum branch length(`minblen`) and/or initial rate omega0.')
		rv[is.na(rv)] <- max(rv, na.rm=T)
		rv <- .hack.times1(rv, td )
	}
	rv
}





# constraints using quad prog 
.optim.Ti5.constrained.limsolve <- function(omegas, td){
		A <- omegas * td$A0 
		B <- td$B0
		B[td$tipEdges] <- td$B0[td$tipEdges] -  unname( omegas[td$tipEdges] * td$sts2 )
		# initial feasible parameter values:
		p0 <- ( coef( lm ( B ~ A -1 , weights = td$W) ) )
		if (any(is.na(p0))){
			warning('Numerical error when performing least squares optimisation. Values are approximate. Try adjusting minimum branch length(`minblen`) and/or initial rate omega0.')
			p0[is.na(p0)] <- max(p0, na.rm=T)
		}
		p1 <- .hack.times1(p0, td )
	
	w <- sqrt(td$W)
	unname( lsei( A = A * w
	 , B = B * w
	 , G = -td$Ain
	 , H = -td$bin
	 , type = 2
	)$X )
}


# <constrained ls>
.optim.Ti2 <- function( omegas, td ){
	success = requireNamespace('mgcv', quietly=TRUE)
	if (!success) stop('The *mgcv* option requires installation of the `mgcv` package. Stopping.') 
		A <- omegas * td$A0 
		B <- td$B0
		B[td$tipEdges] <- td$B0[td$tipEdges] -  unname( omegas[td$tipEdges] * td$sts2 )
		
		# initial feasible parameter values:
		#p0 <- ( coef( lm ( B ~ A -1 , weights = td$W/omegas) ) )
		p0 <- ( coef( lm ( B ~ A -1 , weights = td$W) ) )
		if (any(is.na(p0))){
			warning('Numerical error when performing least squares optimisation. Values are approximate. Try adjusting minimum branch length(`minblen`) and/or initial rate omega0.')
			p0[is.na(p0)] <- max(p0, na.rm=T)
		}
		p1 <- .hack.times1(p0, td )
		
		# design
		M <- list( 
			X  = A
			,p = p1
			,off = c()# rep(0, np)
			,S=list()
			,Ain=-td$Ain
			,bin=-td$bin
			,C=matrix(0,0,0)
			,sp= c()#rep(0,np)
			,y=B
			,w=td$W #/omegas # better performance on lsd tests w/o this  
		)
		o <- mgcv::pcls(M)
	o
}
#</ constrained ls>

.optim.omegas.gammaPoisson1 <- function( Ti, r, gammatheta, td )
{
	blen <- .Ti2blen( Ti, td )
	o <- sapply( 1:nrow(td$tre$edge), function(k){
		sb <- (td$tre$edge.length[k]*td$s)
		lb <- qgamma(1e-6,  shape=r, scale = gammatheta*blen[k] )
		lam_star <- max(lb, gammatheta * blen[k] * (sb + r - 1) / (gammatheta * blen[k] + 1)  )
		ll <- 	dpois( max(0, round(sb)),lam_star, log=T )  + 
		 dgamma(lam_star, shape=r, scale = gammatheta*blen[k], log = T)
		c(lam_star / blen[k] / td$s, ll )
	}) 
	list( omegas = o[1,], ll = unname(sum( o[2,] )), lls  = unname(o[2,])  )
}


.optim.sampleTimes0 <- function( Ti, omegas, estimateSampleTimes, estimateSampleTimes_densities, td, iedge_tiplabel_est_samp_times )
{
	blen <- .Ti2blen( Ti, td )
	o <- sapply( iedge_tiplabel_est_samp_times, function(k) {
		u <- td$tre$edge[k,1]
		v <- td$tre$edge[k,2]
		V <- td$tre$tip.label[v]
		dst <- estimateSampleTimes_densities[[V]]
		tu <- Ti[u-td$n]
		of <- function(tv){
			.blen <- tv - tu 
			-dst(tv, V) -
			  dpois( max(0, round(td$tre$edge.length[k]*td$s))
			  , td$s * .blen * omegas[k] 
			  , log = T)
		}
		lb <- max( tu, estimateSampleTimes[V,'lower'] )
		ub <- max( tu, estimateSampleTimes[V, 'upper'] )
		if (ub == tu & lb == tu) return(tu ) #+ td$minblen
		o  <- optimise(  of, lower = lb, upper = ub )
		o$minimum
	})
	o
}


.mean.rate <- function(Ti, r, gammatheta, omegas, td)
{
	if (is.infinite(r)) return(gammatheta) # poisson model
	blen <- .Ti2blen( Ti, td )
	sum( omegas * blen ) / sum(blen)
}

.hack.times1 <- function(Ti, td)
{
	#~ 	t <- c( td$sts2[td$tre$tip.label], Ti)
	t <- c( td$sts1[td$tre$tip.label], Ti)
	inodes <- (td$n+1):length(t)
	
	repeat
	{
		.t <- t
		.t[inodes] <- pmin(t[inodes], -td$minblen + pmin( t[ td$daughters[inodes,1] ], t[td$daughters[inodes,2] ] ) )
		if (identical( t, .t )) break
		t <- .t
	}
	t[inodes]
}

.dater <- function(tre, sts, s=1e3
 , omega0 = NA
 , minblen = NA
 , maxit=100
 , abstol = .0001
 , searchRoot = 5
 , quiet = TRUE
 , temporalConstraints = TRUE
 , strictClock = FALSE
 , estimateSampleTimes = NULL
 , estimateSampleTimes_densities= list()
 , numStartConditions = 0
 , clsSolver=c('limSolve', 'mgcv')
 , meanRateLimits = NULL
 , ncpu = 1
 , parallel_foreach = FALSE
 , lnd.mean.rate.prior = function(x) 0
 , tiplabel_est_samp_times = NULL
)
{
	clsSolver <- match.arg( clsSolver) 
	# defaults
	CV_LB <- 1e-6 # lsd tests indicate Gamma-Poisson model may be more accurate even in strict clock situation
	cc <- 10
	intree_rooted <- TRUE
	if (!is.rooted(tre)) stop('Tree not rooted.')
	
	EST_SAMP_TIMES  = ifelse( is.null( tiplabel_est_samp_times ),  FALSE, TRUE)
	iedge_tiplabel_est_samp_times <- match( tiplabel_est_samp_times, tre$tip.label[tre$edge[,2]] )
	
	if (is.na(omega0)){
		# guess
		#omega0 <- estimate.mu( tre, sts )
		g0 <- lm(ape::node.depth.edgelength(tre)[1:length(sts)] ~ sts, na.action = na.omit)
		omega0sd <- summary( g0 )$coef[2,2]
		omega0 <- unname( coef(g0)[2] )
		if (omega0 < 0 ){
			warning('Root to tip regression predicts a substition rate less than zero. Tree may be poorly rooted or there may be small temporal signal.')
			omega0 <- abs(omega0)/10
		}
		omega0s <- qnorm( unique(sort(c(.5, seq(.025, .975, l=numStartConditions*2) )))  , omega0, sd = omega0sd )
		omega0s <- omega0s[ omega0s > 0 ]
		#if (!quiet){
		#	cat('initial rates:\n')
		#	print(omega0s)
		#}
	} else{
		omega0s <- c( omega0 )
	}
	if ( any ( omega0s < meanRateLimits[1] ) | any( omega0s > meanRateLimits[2])){
		warning('Initial guess of mean rate falls outside of user-specified limits.')
		omega0s <- omega0s[  omega0s >= meanRateLimits[1] & omega0s <= meanRateLimits[2] ]
	}
	if (length(omega0s)==0){
		warning( 'Setting initial guess of mean rate to be mid-point of *meanRateLimits*')
		omega0s <- (meanRateLimits[1] + meanRateLimits[2]) / 2
	}
	td <- .make.tree.data(tre, sts, s, cc )
	td$minblen <- minblen
	
	omega2ll <- -Inf 
	bestrv <- list()
	for ( omega0 in omega0s ){
		# initial gamma parms with small variance
		r = r0 <- ifelse(strictClock, Inf, sqrt(10))  #sqrt(r) = 10 
		gammatheta = gammatheta0 <- ifelse(strictClock, omega0, omega0 * td$s / r0)
		
		done <- FALSE
		lastll <- -Inf
		iter <- 0
		nEdges <- nrow(tre$edge)
		omegas <- rep( omega0,  nEdges )
		edge_lls <- NA
		rv <- list()
		while(!done){
			if (temporalConstraints){
				if (clsSolver=='limSolve'){
					Ti <- tryCatch( .optim.Ti5.constrained.limsolve ( omegas, td ) 
					 , error = function(e) .optim.Ti2( omegas, td)  )
				} else{
					Ti <- .optim.Ti2( omegas, td)  
				}
			} else{
				Ti <- .optim.Ti0( omegas, td, scale_var_by_rate=FALSE )
			}
			if ( (1 / sqrt(r)) < CV_LB){
				# switch to poisson model
				o <- .optim.omega.poisson0(Ti, .mean.rate(Ti, r, gammatheta, omegas, td), td, lnd.mean.rate.prior , meanRateLimits)
				gammatheta <- unname(o$omega)
				if (!is.infinite(r)) lastll <- -Inf # the first time it switches, do not do likelihood comparison 
				r <- Inf#unname(o$omega)
				ll <- o$ll
				edge_lls <- 0
				omegas <- rep( gammatheta, length(omegas))
			} else{
				o <- .optim.r.gammatheta.nbinom0(  Ti, r, gammatheta, td, lnd.mean.rate.prior )
				r <- o$r
				ll <- o$ll
				gammatheta <- o$gammatheta
				oo <- .optim.omegas.gammaPoisson1( Ti, o$r, o$gammatheta, td ) 
				edge_lls <- oo$lls
				omegas <- oo$omegas
			}
			
			if (EST_SAMP_TIMES)
			{
				o_sts <- .optim.sampleTimes0( Ti, omegas, estimateSampleTimes,estimateSampleTimes_densities, td, iedge_tiplabel_est_samp_times )
				sts[tiplabel_est_samp_times] <- o_sts
				td$sts[tiplabel_est_samp_times] <- o_sts
				td$sts2[tiplabel_est_samp_times] <- o_sts
			}
			
			if (!quiet)
			{ 
				print( data.frame( iteration = iter, median_unadjusted_rate = median(omegas),  coef_of_var =  1 / sqrt(r) ), tmrca = min(Ti), logLik=ll , row.names=iter)
				cat( '---\n' )
			}
			omega <- .mean.rate(Ti, r, gammatheta, omegas, td)
			if ( ll >= lastll ){
				rv <- list( omegas = omegas, r = unname(r), theta = unname(gammatheta), Ti = Ti
				 , meanRate = omega
				 , loglik = ll
				 , edge_lls = edge_lls )
			}
			
			# check convergence
			iter <- iter + 1
			if (iter > maxit) done <- TRUE
			
			if ( abs( ll - lastll ) < abstol) done <- TRUE
			if (ll < lastll) {
				done <- TRUE
			}
			
			lastll <- ll
		}
		if ( omega2ll < rv$loglik){
			bestrv <- rv
			omega2ll <- rv$loglik
		}
	}
	
	rv <- bestrv
	
	#.tre <- tre
	td$minblen <- -Inf; 
	blen <- .Ti2blen( rv$Ti, td )
	#tre$edge.length <- blen 
	#rv$tre <- tre
	
	rv$edge <- tre$edge
	rv$edge.length <- blen
	rv$tip.label <- tre$tip.label
	rv$Nnode <- tre$Nnode
	
	rv$timeOfMRCA <- min(rv$Ti)
	rv$timeToMRCA <- max(sts) - rv$timeOfMRCA
	rv$s <- s
	rv$sts <- sts
	rv$minblen <- minblen
	rv$intree <- tre #.tre
	rv$coef_of_variation <- ifelse( is.numeric(rv$r), 1 / sqrt(rv$r), NA )
	rv$clock <- ifelse( is.infinite(rv$r), 'strict', 'relaxed')
	rv$intree_rooted <- intree_rooted
	rv$is_intree_rooted <- intree_rooted
	rv$temporalConstraints <- temporalConstraints
	rv$estimateSampleTimes <- estimateSampleTimes
	rv$EST_SAMP_TIMES = EST_SAMP_TIMES
	if (!EST_SAMP_TIMES) rv$estimateSampleTimes <- NULL
	rv$estimateSampleTimes_densities <- estimateSampleTimes_densities
	rv$numStartConditions <- numStartConditions
	rv$lnd.mean.rate.prior <- lnd.mean.rate.prior
	rv$meanRateLimits <- meanRateLimits
	
	# add pvals for each edge
	if (rv$clock=='relaxed'){
		rv$edge.p <- with(rv, {
		blen <- pmax(minblen, edge.length)
		ps <- pmin(1 - 1e-12, theta * blen/(1 + theta * blen))
			pnbinom(pmax(0, round(intree$edge.length * s)), size = r, 
				prob = 1 - ps)
		})
	} else{
		rv$edge.p <- with(rv, {
			blen <- pmax(minblen, edge.length)
			ppois(pmax(0, round(intree$edge.length * s)), blen * 
				meanRate * s)
		})
	}
	
	class(rv) <- c('treedater', 'phylo')
	rv
}

#' Estimate a time-scaled tree and fit a molecular clock
#'
#' @details 
#' Estimates the calendar time of nodes in the given phylogenetic
#' tree with branches in units of substitutions per site. The
#' calendar time of each sample must also be specified and the length
#' of the sequences used to estimate the tree. If the tree is not
#' rooted, this function will estimate the root position. 
#' For an introduction to all options and features, see the vignette on Influenza H3N2: vignette("h3n2")
#' 
#' @section References:
#' E.M. Volz and Frost, S.D.W. (2017) Scalable relaxed clock phylogenetic dating. Virus Evolution.
#' 
#' @param tre An ape::phylo which describes the phylogeny with branches in
#'        units of substitutions per site. This may be a rooted or
#'        unrooted tree. If unrooted, the root position will be
#'        estimated by checking multiple candidates chosen by
#'        root-to-tip regression.  If the tree has multifurcations,
#'        these will be resolved and a binary tree will be returned.
#' @param sts Vector of sample times for each tip in phylogenetic tree.
#'        Vector must be named with names corresponding to
#'        tre$tip.label.
#' @param s Sequence length (numeric). This should correspond to sequence length used in phylogenetic analysis and will not necessarily be the same as genome length. 
#' @param omega0 Initial guess of the mean substitution rate (substitutions
#'        per site per unit time). If not provided, will guess using
#'        root to tip regression.
#' @param minblen Minimum branch length in calendar time. By default, this will
#'        be the range of sample times (max - min) divided by sample
#'        size.
#' @param maxit Maximum number of iterations 
#' @param abstol Difference in log likelihood between successive iterations for convergence.
#' @param searchRoot Will search for the optimal root position using the top
#'        matches from root-to-tip regression.  If searchRoot=x, dates
#'        will be estimated for x trees, and the estimate with the
#'        highest likelihood will be returned.
#' @param quiet If TRUE, will suppress messages during execution 
#' @param temporalConstraints  If TRUE, will enforce the condition that an
#'        ancestor node in the phylogeny occurs before all progeny.
#'        Equivalently, this will preclude negative branch lengths.
#'        Note that execution is faster if this option is FALSE.
#' @param strictClock If TRUE, will fit a Poisson evolutionary model without
#'        rate variation.
#' @param estimateSampleTimes If some sample times are not known with certainty,
#'         bounds can be provided with this option. This should take the
#'         form of a data frame with columns 'lower' and 'upper'
#'         providing the sample time bounds for each uncertain tip. Row
#'         names of the data frame should correspond to elements in
#'         tip.label of the input tree. Tips with sample time bounds in
#'         this data frame do not need to appear in the *sts* argument,
#'         however if they are included in *sts*, that value will be
#'         used as a starting condition for optimisation.
#' @param estimateSampleTimes_densities An optional named list of log densities
#'           which would be used as priors for unknown sample times. Names
#'           should correspond to elements in tip.label with uncertain
#'           sample times.
#' @param numStartConditions Will attempt optimisation from more than one starting point if >0
#' @param clsSolver Which package should be used for constrained least-squares? Options are "mgcv" or "limSolve"
#' @param meanRateLimits Optional constraints for the mean substitution rate 
#' @param ncpu Number of threads for parallel computing 
#' @param parallel_foreach If TRUE, will use the "foreach" package instead of the "parallel" package. This may work better on some HPC systems. 
#' 
#' @return A time-scaled tree and estimated molecular clock rate 
#' 
#' @author Erik M Volz <erik.volz@gmail.com>
#' 
#' @seealso 
#' ape::chronos
#' ape::estimate.mu
#'
#' @examples
#' ## simulate a random tree and sample times for demonstration
#' # make a random tree:
#' tre <- ape::rtree(50)
#' # sample times based on distance from root to tip:
#' sts <- setNames(  
#'   ape::dist.nodes( tre)[(length(tre$tip.label)+1), 1:(length(tre$tip.label)+1)]
#'   , tre$tip.label)
#' # modify edge length to represent evolutionary distance with rate 1e-3:
#' tre$edge.length <- tre$edge.length * 1e-3
#' # treedater: 
#' td <- dater( tre, sts =sts )
#'
#'
#' @export 
dater <- function(tre, sts, s=1e3
 , omega0 = NA
 , minblen = NA
 , maxit=100
 , abstol = .0001
 , searchRoot = 5
 , quiet = TRUE
 , temporalConstraints = TRUE
 , strictClock = FALSE
 , estimateSampleTimes = NULL
 , estimateSampleTimes_densities= list()
 , numStartConditions = 0
 , clsSolver=c('limSolve', 'mgcv')
 , meanRateLimits = NULL
 , ncpu = 1
 , parallel_foreach = FALSE
)
{ 
	clsSolver <- match.arg( clsSolver) #clsSolver[1]	
	if (!is.binary( tre ) ){
		cat( 'Note: *dater* called with non binary tree. Will proceed after resolving polytomies.\n' )
		if ( !is.rooted( tre )){
			tre <- unroot( multi2di( tre ) ) 
		} else{
			tre <-  multi2di( tre ) 
		}
	}
	
	if (class(tre)[1]=='treedater'){
		cat('Note: *dater* called with treedater input tree. Will use rooted tree with branch lengths in substitions.\n')
		tre <- tre$intree
	}
	
	# optional limits or prior for mean rate
	lnd.mean.rate.prior <- function(x) 0 #dunif( x , 0, Inf, log=TRUE )
	if (is.null( meanRateLimits) ) {
		meanRateLimits <- c( 0, Inf)
	} else if ( class(meanRateLimits)=='function'){
		lnd.mean.rate.prior <- meanRateLimits
		meanRateLimits <- c(0, Inf)
	} else{
		MEANRATEERR <- '*meanRateLimits* should be a length 2 vector providing bounds on the mean rate parameter OR a function providing the log prior density of the mean rate parameter. '
		if (length( meanRateLimits) != 2) stop(MEANRATEERR)
		if (meanRateLimits[2] <= meanRateLimits[1]) stop(MEANRATEERR)
		#lnd.mean.rate.prior <- function(x) dunif( x , meanRateLimits[1], meanRateLimits[2], log= TRUE )
		lnd.mean.rate.prior <- function(x) ifelse( x >= meanRateLimits[1] & x <= meanRateLimits[2], 0, -Inf )
	}
	
	numStartConditions <- max(0, round( numStartConditions )) # number of omega0 to try for optimisation
	
	EST_SAMP_TIMES <- TRUE
	EST_SAMP_TIMES_ERR <- 'estimateSampleTimes must specify a data frame with tip.label as row names and with columns `upper` and `lower`. You may also provide a named list of log density functions (improper priors for sample times).\n'
	if (is.null(estimateSampleTimes)) EST_SAMP_TIMES <- FALSE
	tiplabel_est_samp_times <- NULL 
	if (EST_SAMP_TIMES){
		if (class(estimateSampleTimes)=='data.frame'){
			if ( !('lower' %in% colnames(estimateSampleTimes)) | !('upper' %in% colnames(estimateSampleTimes) ) ){
				stop(EST_SAMP_TIMES_ERR)
			}
			if ( any (estimateSampleTimes$lower > estimateSampleTimes$upper) ){
				stop(EST_SAMP_TIMES_ERR )
			}
			estimateSampleTimes <- estimateSampleTimes[ estimateSampleTimes$lower < estimateSampleTimes$upper ,]
			for (tl in rownames(estimateSampleTimes)){
				if (!(tl %in% names( estimateSampleTimes_densities))){
					estimateSampleTimes_densities[[tl]] <-  function(x,tl) dunif(x, min= estimateSampleTimes[tl,'lower'], max=estimateSampleTimes[tl,'upper'] , log = TRUE) #
				}
			}
		} else {
			stop(EST_SAMP_TIMES_ERR)
		}
		tiplabel_est_samp_times <- intersect( rownames(estimateSampleTimes), tre$tip.label)
		iedge_tiplabel_est_samp_times <- match( tiplabel_est_samp_times, tre$tip.label[tre$edge[,2]] )
	}
		
	# check for missing sample times, impute missing if needed 
	#if (any(is.na(sts))) stop( 'Some sample times are NA.' )
	stinfo_provided <- union( names(na.omit(sts)), rownames(estimateSampleTimes))
	stinfo_not_provided <-   setdiff( tre$tip.label, stinfo_provided ) 
	if (length( stinfo_not_provided ) > 0){
		cat( 'NOTE: Neither sample times nor sample time bounds were provided for the following lineages:\n')
		cat( stinfo_not_provided )
		cat('\n Provide sampling info or remove these lineages from the tree. Stopping.\n ') 
		stop('Missing sample time information.' )
	}
	initial_st_should_impute <- setdiff( rownames(estimateSampleTimes), names(na.omit(sts)))
	if (length( initial_st_should_impute ) > 0){
		cat('NOTE: initial guess of sample times for following lineages was not provided:\n')
		cat ( initial_st_should_impute ) 
		cat('\n') 
		cat( 'Will proceed with midpoint of provided range as initial guess of these sample times.\n')
		sts[initial_st_should_impute] <- rowMeans( estimateSampleTimes )[initial_st_should_impute] 
	}
	
	if (is.null(names(sts))){
		if (length(sts)!=length(tre$tip.label)) stop('Sample time vector length does not match number of lineages.')
		names(sts) <- tre$tip.label
	}
	sts <- sts[tre$tip.label]
	
	if (is.na(minblen)){
		minblen <- diff(range(sts))/ length(sts) / 10 #TODO choice of this parm is difficult, may require sep optim / crossval
		cat(paste0('Note: Minimum temporal branch length  (*minblen*) set to ', minblen, '. Increase *minblen* in the event of convergence failures. \n'))
	}
	if (!is.na(omega0) & numStartConditions > 0 ){
		warning('omega0 provided incompatible with numStartConditions > 0. Setting numStartConditions to zero.')
		numStartConditions <- 0
	}
	
	intree_rooted <- TRUE
	if (!is.rooted(tre)){
		intree_rooted <- FALSE
		cat( 'Tree is not rooted. Searching for best root position. Increase searchRoot to try harder.\n')
		searchRoot <- round( searchRoot )
		rtres <- .multi.rtt(tre, sts, topx=searchRoot, ncpu = ncpu)
		if (ncpu > 1 )
		{
			if (parallel_foreach){
				success = requireNamespace('foreach', quietly=TRUE)
				if (!success) stop('*parallel_foreach* requires the `foreach` package. Stopping.')
				`%dopar%` <- foreach::`%dopar%`
				tds <- foreach::foreach( t = iterators::iter( rtres )) %dopar% {
					.dater( t, sts, s = s, omega0=omega0, minblen=minblen, maxit=maxit,abstol=abstol
						, strictClock = strictClock, temporalConstraints = temporalConstraints, quiet = quiet
						, estimateSampleTimes = estimateSampleTimes
						, estimateSampleTimes_densities = estimateSampleTimes_densities  
						, numStartConditions = numStartConditions
						, meanRateLimits = meanRateLimits
						, lnd.mean.rate.prior =  lnd.mean.rate.prior 
						, tiplabel_est_samp_times = tiplabel_est_samp_times
						) 
				}
			} else{
				tds <- parallel::mclapply( rtres, function(t){
					.dater( t, sts, s = s, omega0=omega0, minblen=minblen, maxit=maxit,abstol=abstol
						, strictClock = strictClock, temporalConstraints = temporalConstraints, quiet = quiet
						, estimateSampleTimes = estimateSampleTimes
						, estimateSampleTimes_densities = estimateSampleTimes_densities  
						, numStartConditions = numStartConditions
						, meanRateLimits = meanRateLimits
						, lnd.mean.rate.prior =  lnd.mean.rate.prior 
						, tiplabel_est_samp_times = tiplabel_est_samp_times
						) 
				}, mc.cores = ncpu )
			}
		} else{
			tds <- lapply( rtres, function(t) {
				.dater( t, sts, s = s, omega0=omega0, minblen=minblen, maxit=maxit,abstol=abstol
					, strictClock = strictClock, temporalConstraints = temporalConstraints, quiet = quiet
					, estimateSampleTimes = estimateSampleTimes
					, estimateSampleTimes_densities = estimateSampleTimes_densities  
					, numStartConditions = numStartConditions
					, meanRateLimits = meanRateLimits
					, lnd.mean.rate.prior = lnd.mean.rate.prior 
					, tiplabel_est_samp_times = tiplabel_est_samp_times
					) 
			})
		}
		lls <- sapply( tds, function(td) td$loglik )
		td <- tds [[ which.max( lls ) ]]
		td$intree_rooted <- FALSE
		return ( td )
	} else{
		cat( 'Tree is rooted. Not estimating root position.\n')
	}
	.dater( tre, sts, s = s, omega0=omega0, minblen=minblen, maxit=maxit,abstol=abstol
		, strictClock = strictClock, temporalConstraints = temporalConstraints
		, quiet = quiet
		, estimateSampleTimes = estimateSampleTimes
		, estimateSampleTimes_densities = estimateSampleTimes_densities  
		, numStartConditions = numStartConditions
		, meanRateLimits = meanRateLimits
		, lnd.mean.rate.prior = lnd.mean.rate.prior 
		, tiplabel_est_samp_times = tiplabel_est_samp_times
	)
}

#' @export 
print.treedater <- function(x, ...){
    cl <- oldClass(x)
    oldClass(x) <- cl[cl != "treedater"]
    print(x$intree)
    cat('\n Time of common ancestor \n' )
    cat(paste( x$timeOfMRCA, '\n') )
    cat('\n Time to common ancestor (before most recent sample) \n' )
    cat(paste( x$timeToMRCA, '\n') )
    cat( '\n Mean substitution rate \n')
    cat(paste( x$meanRate , '\n'))
    cat( '\n Strict or relaxed clock \n')
    cat(paste( x$clock , '\n'))
    cat( '\n Coefficient of variation of rates \n')
    cat(paste( x$coef_of_variation, '\n' ))
    
    invisible(x)
}

#' @export 
summary.treedater <- function(object, ...) {
    stopifnot(inherits(object, "treedater"))
    print.treedater( object )
}

#' Produce a goodness of fit plot
#'
#' The sorted tail probabilties (p values) for each edge in the tree under the fitted model
#'
#' @param td A treedater object generated by the \code{dater} function
#'
#' @export 
goodnessOfFitPlot <- function(td)
{
	stopifnot(inherits(td, "treedater"))
with( td, 
	{
		plot( 1:length(edge.p)/length(edge.p), sort (edge.p ) , type = 'l', xlab='Theoretical quantiles', ylab='Edge p value'); 
		abline( a = 0, b = 1 )
	})
}

