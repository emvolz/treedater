#~ Treedater: fast relaxed molecular clock dating
#~     Copyright (C) 2025  Erik Volz
#~     This program is free software: you can redistribute it and/or modify
#~     it under the terms of the GNU General Public License as published by
#~     the Free Software Foundation, either version 3 of the License, or
#~     (at your option) any later version.


#' Compute a vector of numeric sample times from labels in a sequence alignment or phylogeny
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

	postorder_internal_nodes <- unique( tre$edge[ ape::postorder( tre ),1] )

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
	 , postorder_internal_nodes = postorder_internal_nodes
	)
}


.optim.r.gammatheta.nbinom0 <- function(  Ti, r0, gammatheta0, td, lnd.mean.rate.prior)
{
	blen <- .Ti2blen( Ti, td )
	MINR <- 1e-3
	MAXR <- 1e3 
	#NOTE relative to wikipedia page on NB:
	#size = r
	#1-prob = p
	of <- function(x)
	{
		r <- min(MAXR, max(MINR, exp( x['lnr'] ) ) )
		gammatheta <- exp( x['lngammatheta'] )
		if (is.infinite(r) | is.infinite(gammatheta)) return(Inf)
		ps <- pmin(1 - 1e-5, gammatheta*blen / ( 1+ gammatheta * blen ) )
		ov <- -sum( dnbinom( pmax(0, round(td$tre$edge.length*td$s))
		  , size= r, prob=1-ps,  log = TRUE) )
		mr <-  r * gammatheta / td$s
		ov <- ov - unname(lnd.mean.rate.prior( mr ))
		ov
	}
	x0 <- c( lnr = unname(log(r0)), lngammatheta = unname( log(gammatheta0)))
	#o <- optim( par = x0, fn = of, method = 'BFGS' )
	if ( is.infinite( of( x0 )) ) stop( 'Can not optimize rate parameters from initial conditions. Try adjusting *meanRateLimits*.' )
	o <- optim( par = x0, fn = of)
	r <- min(MAXR, max(MINR, unname( exp( o$par['lnr'] )) ))
	gammatheta <- unname( exp( o$par['lngammatheta'] ))
	list( r = r, gammatheta=gammatheta, ll = -o$value)
}

# additive varianc model of x didelot
.optim.nbinom1 <- function(  Ti, mu0, sp0, td, lnd.mean.rate.prior)
{
	blen <- .Ti2blen( Ti, td )

	#NOTE relative to wikipedia page on NB:
	#size = r
	#1-prob = p
	of <- function(x)
	{
		sp <- exp( x[ 'lnsp' ] )
		mu <- exp( x['lnmu'] )
		if (is.infinite(sp) | is.infinite(mu)) return(Inf)
		sizes <- mu * blen / sp
		ov <- -sum( dnbinom( pmax(0, round(td$tre$edge.length*td$s))
		  , size = sizes
		  , prob = 1 - sp / ( 1+sp ),  log = TRUE) )
		ov <- ov - unname(lnd.mean.rate.prior( mu  / td$s))
		ov
	}
	x0 <- c( lnmu = unname(log(mu0)), lnsp = unname( log(sp0)))
	#o <- optim( par = x0, fn = of, method = 'BFGS' )
	if ( is.infinite( of( x0 )) ) stop( 'Can not optimize rate parameters from initial conditions. Try adjusting *meanRateLimits*.' )
	o <- optim( par = x0, fn = of)
	mu <- unname( exp( o$par['lnmu'] ))
	sp <- unname( exp( o$par['lnsp'] ))
	list( mu = mu, sp = sp, ll = -o$value)
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


# using additive relaxed clock model
.optim.omegas.gammaPoisson2 <- function( Ti, mu, sp , td )
{
	blen <- .Ti2blen( Ti, td )
	o <- sapply( 1:nrow( td$tre$edge ), function(k) {
		tau <- blen[k]
		x <- (td$tre$edge.length[k]*td$s)
		#lb <- qgamma(1e-6,  shape= mu*tau/sp , scale = sp ) # gives values too close to zero
		lb <- mu / 100# TODO would be nice to have a nice way to  choose this lower bound
		lam_star <- max(lb
		 , (mu*tau - sp + x * sp)  /  (sp + 1 )  # yup.
		# , (mu*tau - sp + x * sp)  /  (sp + tau ) # nope.
		)
		ll <- dpois( max(0, round(x)),lam_star, log=T )  +
		 dgamma(lam_star, shape=mu*tau/sp , scale = sp / tau , log = T)
		c( lam_star / tau / td$s, ll )
	})
	list( omegas = o[1,] , ll = unname( sum( o[2,] )), lls = unname(o[2,] ) )
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

.Ti2blen <- function(Ti, td, enforce_minblen = TRUE ){
	Ti <- c( td$sts[td$tre$tip.label], Ti)
	elmat <- cbind( Ti[td$tre$edge[,1]], Ti[td$tre$edge[,2]])
	if ( enforce_minblen )
		return( pmax(td$minblen,-elmat[,1]  + elmat[,2] ) )
	else
		return( -elmat[,1]  + elmat[,2] )
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


.fix.Ti <- function(Ti, td){
	.Ti <- Ti
	n <- 1 + length(Ti)
	sts_Ti <- c( td$sts1[td$tre$tip.label], Ti)
	for ( i in td$postorder_internal_nodes-n ){
		uv <- td$daughters[ i+n, ]
		if ( !is.na(sts_Ti[uv[1]]))
			Ti[i] <- min( Ti[i], sts_Ti[uv[1]] - td$minblen  )
		if ( !is.na(sts_Ti[uv[2]]))
			Ti[i] <- min( Ti[i], sts_Ti[uv[2]] - td$minblen  )
		sts_Ti[ i + n ] <- Ti[i]
	}
	Ti
}


.dater <- function(tre, sts, s=1e3
 , omega0 = NA
 , minblen = NA
 , maxit=100
 , abstol = .0001
 , searchRoot = 5
 , quiet = TRUE
 , temporalConstraints = TRUE
 , clock = c('strict',  'uncorrelated', 'additive')
 , estimateSampleTimes = NULL
 , estimateSampleTimes_densities= list()
 , numStartConditions = 1
 , clsSolver=c('limSolve', 'mgcv')
 , meanRateLimits = NULL
 , ncpu = 1
 , parallel_foreach = FALSE
 , lnd.mean.rate.prior = function(x) 0
 , tiplabel_est_samp_times = NULL
)
{
	clsSolver <- match.arg( clsSolver, choices = c('limSolve', 'mgcv'))
	clock <- match.arg( clock , choices = c('uncorrelated', 'additive', 'strict') )
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
		# start conditiosn based on rtt
			g0 <- lm(ape::node.depth.edgelength(tre)[1:length(sts)] ~ sts, na.action = na.omit)
			omega0sd <- suppressWarnings(summary( g0 ))$coef[2,2]
			omega0 <- unname( coef(g0)[2] )
			if (omega0 < 0 ){
				warning('Root to tip regression predicts a substition rate less than zero. Tree may be poorly rooted or there may be small temporal signal.')
				omega0 <- abs(omega0)/10
			}
			omega0s <- qnorm( unique(sort(c(.5, seq(.001, .999, l=numStartConditions*2) )))  , omega0, sd = omega0sd )
		#start conditiosn based on earliest sample
			D <- ape::cophenetic.phylo( tre ) [1:ape::Ntip(tre), 1:ape::Ntip(tre)]
			esi <- which.min( sts )[1]
			g1 <- lm( D[esi,] ~ sts )
			omega0sd.1 <- suppressWarnings(summary( g1 ))$coef[2,2]
			omega0.1 <- unname( coef(g1)[2] )
			omega0s.1 <- qnorm( unique(sort(c(.5, seq(.001, .999, l=numStartConditions*2) )))  , omega0.1, sd = omega0sd.1 )


		omega0s <- sort( unique( c( omega0s, omega0s.1 ) ))
		omega0s <- omega0s[ omega0s > 0 ]
		cat( paste( 'Initial guesses of substitution rate:', paste(collapse=',', omega0s), '\n')  )
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
		r = r0 <- ifelse(clock=='strict', Inf, sqrt(10))  #sqrt(r) = 10
		gammatheta = gammatheta0 <- ifelse(clock=='strict', omega0, omega0 * td$s / r0)
		mu = omega0  * td$s
		sp = 1e-2

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
			Ti <- .fix.Ti( Ti, td )
			if ( (1 / sqrt(r)) < CV_LB){
				# switch to poisson model
				o <- .optim.omega.poisson0(Ti
				  , mean(omegas)
				  , td, lnd.mean.rate.prior , meanRateLimits)
				gammatheta <- unname(o$omega)
				if (!is.infinite(r)) lastll <- -Inf # the first time it switches, do not do likelihood comparison
				r <- Inf#unname(o$omega)
				ll <- o$ll
				edge_lls <- 0
				omegas <- rep( gammatheta, length(omegas))
			} else{
				if (clock == 'uncorrelated')
				{
					o <- .optim.r.gammatheta.nbinom0(  Ti, r, gammatheta, td, lnd.mean.rate.prior )
					r <- o$r
					ll <- o$ll
					gammatheta <- o$gammatheta
					oo <- .optim.omegas.gammaPoisson1( Ti, o$r, o$gammatheta, td )
				} else if (clock=='additive'){
					o <- .optim.nbinom1( Ti, mu, sp, td, lnd.mean.rate.prior )
					mu = o$mu
					sp = o$sp
					ll = o$ll
					oo = .optim.omegas.gammaPoisson2( Ti, mu, sp , td )
				} else{
					stop('invalid value for *clock*')
				}
				edge_lls <- oo$lls
				omegas <- oo$omegas
			}

			if (EST_SAMP_TIMES)
			{
				o_sts <- .optim.sampleTimes0( Ti, omegas, estimateSampleTimes,estimateSampleTimes_densities, td, iedge_tiplabel_est_samp_times )
				sts[tiplabel_est_samp_times] <- o_sts
				td$sts[tiplabel_est_samp_times] <- o_sts
				td$sts1[tiplabel_est_samp_times] <- o_sts
				td$sts2[tiplabel_est_samp_times] <- o_sts
			}

			if (!quiet)
			{
				print( data.frame( iteration = iter, median_unadjusted_rate = median(omegas)
				 ,  coef_of_var =  sd(omegas) / mean(omegas) # 1 / sqrt(r)
				 , tmrca = min(Ti), logLik=ll , row.names=iter) )
				cat( '---\n' )
			}

			if (clock !='strict'){
				blen <- .Ti2blen( Ti, td )
				omega <- sum( omegas * blen ) / sum(blen)
			} else{#strict
				omega <- omegas[1]
			}

			if ( clock=='strict'){
				meanRate = omegas[1]
			} else if (clock=='uncorrelated' ){
				meanRate <- ( r * gammatheta / td$s )
				#NB:  r, tau phi / (1 = tau * phi )
				# Gamma: r ,  phi * tau
			} else if( clock=='additive'){
				meanRate <- mu / td$s
			}


			if ( ll >= lastll ){
				rv <- list( omegas = omegas, r = unname(r), theta = unname(gammatheta), Ti = Ti
				 , adjusted.mean.rate = omega
				 , mean.rate = meanRate
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

	# one last round of st estimation b/c td$sts may not match bestrv
	if (EST_SAMP_TIMES)
	{
		o_sts <- .optim.sampleTimes0( rv$Ti, rv$omegas, estimateSampleTimes,estimateSampleTimes_densities, td, iedge_tiplabel_est_samp_times )
		sts[tiplabel_est_samp_times] <- o_sts
		td$sts[tiplabel_est_samp_times] <- o_sts
		td$sts1[tiplabel_est_samp_times] <- o_sts
		td$sts2[tiplabel_est_samp_times] <- o_sts
	}

	blen <- .Ti2blen( rv$Ti, td , enforce_minblen=FALSE)
	rv$edge <- tre$edge
	rv$edge.length <- blen
	rv$tip.label <- tre$tip.label
	rv$Nnode <- tre$Nnode

	#check whether node.label exists. If node.label exists in a tree
	#then it will be copied to rv
	if(!is.null(tre$node.label)){
	  rv$node.label <- tre$node.label
	}

	rv$timeOfMRCA <- min(rv$Ti)
	rv$timeToMRCA <- max(sts) - rv$timeOfMRCA
	rv$s <- s
	rv$sts <- sts
	rv$minblen <- minblen
	rv$intree <- tre #.tre
	rv$coef_of_variation <- sd(rv$omegas) / mean(rv$omegas) #ifelse( is.numeric(rv$r), 1 / sqrt(rv$r), NA )  #
	rv$clock <- clock
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
	rv$relaxedClockModel = clock
	rv$mu <- mu
	rv$sp <- sp
	rv$omega0 = omega0
	rv$omega0s = omega0s


	# add pvals for each edge
	if (rv$clock=='uncorrelated'){
		rv$edge.p <- with(rv, {
			blen <- pmax(minblen, edge.length)
			ps <- pmax(1e-12, pmin(1 - 1e-12, theta * blen/(1 + theta * blen)) )
			pnbinom(pmax(0, round(intree$edge.length * s)), size = r,
					prob = 1 - ps)
			})
	} else if ( rv$clock=='additive'){
			rv$edge.p <- with(rv, {
				blen <- pmax(minblen, edge.length)
				sizes = mu * blen / sp
				pnbinom(pmax(0, round(intree$edge.length * s)), size = sizes,
					prob = 1 - sp / (1+sp) )
			})
	} else if (rv$clock=='strict') {
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
#' Multiple molecular clock models are supported including a strict clock and two variations on relaxed clocks. The 'uncorrelated' relaxed clock is the Gamma-Poisson mixture presented by Volz and Frost (2017), while the 'additive' variance model was developed by Didelot & Volz (2019).
#'
#' @section References:
#' Volz, E. M., & Frost, S. D. (2017). Scalable relaxed clock phylogenetic dating. Virus evolution, 3(2), vex025.
#' 
#' Didelot, X., Siveroni, I., & Volz, E. M. (2021). Additive uncorrelated relaxed clock models for the dating of genomic epidemiology phylogenies. Molecular Biology and Evolution, 38(1), 307-317.
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
#' @param omega0 Vector providing initial guess or guesses of the mean substitution rate (substitutions
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
#' @param clock The choice of molecular clock model. Choices are 'strict'(default), 'uncorrelated', or 'additive'
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
#' sts <- setNames( ape::node.depth.edgelength( tre )[1:ape::Ntip(tre)], tre$tip.label)
#' # modify edge length to represent evolutionary distance with rate 1e-3:
#' tre$edge.length <- tre$edge.length * 1e-3
#' # treedater:
#' td <- dater( tre, sts =sts , s = 1000, clock='strict', omega0=.0015)
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
 , clock = c('strict' , 'uncorrelated', 'additive')
 , estimateSampleTimes = NULL
 , estimateSampleTimes_densities= list()
 , numStartConditions = 1
 , clsSolver=c('limSolve', 'mgcv')
 , meanRateLimits = NULL
 , ncpu = 1
 , parallel_foreach = FALSE
)
{
	clsSolver <- match.arg( clsSolver, choices = c('limSolve', 'mgcv'))
	clock <- match.arg( clock , choices = c('strict' , 'uncorrelated', 'additive') )
	# defaults
	if ( is.na( omega0 ) ){
		cat('Note: Initial guess of substitution rate not provided. Will attempt to guess starting conditions. Provide initial guesses of the rate using *omega0* parameter. \n')
	}
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
	} else if ( is.function(meanRateLimits)){
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
		if (is.data.frame(estimateSampleTimes)){
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
						, clock = clock
						, temporalConstraints = temporalConstraints, quiet = quiet
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
						, clock = clock
						, temporalConstraints = temporalConstraints, quiet = quiet
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
					, clock = clock
					, temporalConstraints = temporalConstraints, quiet = quiet
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
		.fitDiagnostics ( td )
		return ( td )
	} else{
		cat( 'Tree is rooted. Not estimating root position.\n')
	}
	td = .dater( tre, sts, s = s, omega0=omega0, minblen=minblen, maxit=maxit,abstol=abstol
		, clock = clock
		, quiet = quiet
		, estimateSampleTimes = estimateSampleTimes
		, estimateSampleTimes_densities = estimateSampleTimes_densities
		, numStartConditions = numStartConditions
		, meanRateLimits = meanRateLimits
		, lnd.mean.rate.prior = lnd.mean.rate.prior
		, tiplabel_est_samp_times = tiplabel_est_samp_times
	)
	.fitDiagnostics ( td )
	return( td )
}

#' @export
print.treedater <- function(x, ...){
	.fitDiagnostics( x )

    cl <- oldClass(x)
    oldClass(x) <- cl[cl != "treedater"]
    print(x$intree)
    cat('\n Time of common ancestor \n' )
    cat(paste( x$timeOfMRCA, '\n') )
    cat('\n Time to common ancestor (before most recent sample) \n' )
    cat(paste( x$timeToMRCA, '\n') )
    cat( '\n Weighted mean substitution rate (adjusted by branch lengths) \n')
    cat(paste( x$adjusted.mean.rate , '\n'))
    cat( '\n Unadjusted mean substitution rate \n')
    cat(paste( x$mean.rate , '\n'))
    cat( '\n Clock model  \n')
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
#' The sorted tail probabilities (p values) for each edge in the tree under the fitted model
#'
#' @param td A treedater object generated by the \code{dater} function
#'
#' @export
goodnessOfFitPlot <- function(td)
{
	stopifnot(inherits(td, "treedater"))
with( td,
	{
		graphics::plot( 1:length(edge.p)/length(edge.p), sort (edge.p ) , type = 'l', xlab='Theoretical quantiles', ylab='Edge p value', xlim = c(0,1), ylim = c(0,1));
		abline( a = 0, b = 1 )
	})
}

# help message when model fit is poor
.TROUBLESHOOT0 <- '
The following steps may help to fix or alleviate common problems:
* Check that the vector of sample times is correctly named and that the units are correct.
* If passing a rooted tree, make sure that the root position was chosen correctly, or estimate the root position by passing an unrooted tree (e.g. pass ape::unroot(tree))
* The root position may be poorly estimated. Try increasing the _searchRoot_ parameter in order to test more lineages as potential root positions.
* The model may be fitted by a relaxed or strict molecular clock. Try changing the _clock_ parameter
* A poor fit may be due to a small number of lineages with unusual / outlying branch lengths which can occur due to sequencing error or poor alignment. Try the *outlierTips* command to identify and remove these lineages.
* Check that there is adequate variance in sample times in order to estimate a molecular clock by doing a root-to-tip regression. Try the *rootToTipRegressionPlot* command. If the clock rate can not be reliably estimated, you can fix the value to a range using the _meanRateLimits_ option which would estimate a time tree given the previous estimate of clock rates.
'

# Check for common problems in treedater fit and suggest solutions if applicable
.fitDiagnostics <- function( td ){
	stopifnot(inherits(td, "treedater"))
	p <- td$edge.p
	cv <- td$coef_of_variation
	cvprob = FALSE
	pprob = FALSE
	if (cv > 1 ){
		cat('
NOTE: The estimated coefficient of variation of clock rates is high (>1). Sometimes this indicates a poor fit and/or a problem with the data.
		\n')
		cvprob <- TRUE
	}
	#~ success = requireNamespace('harmonicmeanp', quietly=TRUE)
	#~ 	if ( success ){
	#~ 	  pp <- harmonicmeanp::p.hmp( p )
	#~ 	} else {
	#~ 	  pp <- min( p.adjust( p, 'BH' ))
	#~ 	}
	pp <- min( p.adjust( p, 'BH' ))

	if ( pp < .025 ){
		pprob <-  TRUE
		cat( '
NOTE: The p values for lineage clock rates show at least one outlying value after adjusting for multiple-testing.  This indicates a poor fit to the data for at least a portion of the phylogenetic tree. To visualize the distribution p-values, use *goodnessOfFitPlot*.
		\n')
	}
	if (pprob | cvprob ){
		cat( .TROUBLESHOOT0 )
	}
}


#' Plot evolutionary distance from root to sample times and estimated internal node times and regression lines
#'
#' If a range of sample times was given, these will be estimated. Red and black respectively indicate sample and internal nodes.
#' This function will print statistics computed from the linear regression model.
#'
#' @param td A fitted treedater object
#' @param show.tip.labels If TRUE, the names of each sample will be plotted at the their corresponding time and evolutionary distance
#' @param textopts An optional list of parameters for plotted tip labels. Passed to the *text* function.
#' @param pointopts An optional list of parameters for plotted points if showing tip labels. Passed to the *points* function.
#' @param ... Additional arguments are passed to plot
#' @return The fitted linear model (class 'lm')
#' @examples
#' ## simulate a random tree and sample times for demonstration
#' # make a random tree:
#' tre <- ape::rtree(50)
#' # sample times based on distance from root to tip:
#' sts <- setNames( ape::node.depth.edgelength( tre )[1:ape::Ntip(tre)], tre$tip.label)
#' # modify edge length to represent evolutionary distance with rate 1e-3:
#' tre$edge.length <- tre$edge.length * 1e-3
#' # treedater:
#' td <- dater( tre, sts =sts, clock='strict', s = 1000, omega0=.0015 )
#' # root to tip regression:
#' fit = rootToTipRegressionPlot( td )
#' summary(fit)
#'
#' @export
rootToTipRegressionPlot <- function(td, show.tip.labels=FALSE, textopts = NULL, pointopts=NULL, ... ){
	stopifnot( inherits( td, 'treedater'))
	dT <- ape::node.depth.edgelength( td  )
	dG <- ape::node.depth.edgelength( td$intree )
	#scatter.smooth( dT, dG )
	sts <- (td$timeOfMRCA+dT[1:ape::Ntip(td)])
	nts <- (td$timeOfMRCA+dT)
	mtip  <- lm( dG[1:ape::Ntip(td)] ~ sts )
	mall  <- lm( dG ~ nts )
	if ( !show.tip.labels){
		graphics::plot( dT + td$timeOfMRCA, dG
		  , col = c(rep('red', ape::Ntip(td)), rep('black', ape::Nnode(td) ) )
		  , xlab = ''
		  , ylab = 'Evolutionary distance'
		  , ...
		)
	} else {
		i <- 1:ape::Ntip(td)
		j <- (ape::Ntip(td)+1):( ape::Ntip(td) + ape::Nnode(td) )
		graphics::plot( x = NULL, y = NULL
		  , xlab = ''
		  , ylab = 'Evolutionary distance'
		  , xlim = range( dT + td$timeOfMRCA)
		  , ylim = range( dG )
		  , ...
		)
		do.call( graphics::points, c( pointopts, list(x =   dT[j] + td$timeOfMRCA, y = dG[j] )))
		do.call( graphics::text, c( list( x = dT[i] + td$timeOfMRCA, y = dG[i] , labels=td$tip.label ), textopts ) )
	}
	graphics::abline( a = coef(mtip)[1], b = coef(mtip)[2], col = 'red' )
	graphics::abline( a = coef(mall)[1], b = coef(mall)[2], col = 'black' )
	smtip <- summary( mtip )
	cat(paste( 'Root-to-tip mean rate:', coef(mtip)[2], '\n'))
	cat(paste( 'Root-to-tip p value:', smtip$coefficients[2, 4 ], '\n'))
	cat(paste( 'Root-to-tip R squared (variance explained):', smtip$r.squared, '\n'))
	cat('Returning fitted linear model.\n')
	invisible(mtip)
}


#' Sample node dates conditional on root time using Gibbs sampling and date prior
#'
#' This function is useful for 'smoothing out' time trees that have many adjacent small branch lengths (essential polytomies). It returns a list of smoothed trees.
#'
#' @param dtr A treedater fit
#' @param iter Number of iterations (every node time is resampled in each iteration )
#' @param burn_pc Remove this proportion as burn-in
#' @param returnTrees Integer number of trees to return.
#' @param res Time resolution for proposing new node times
#' @param report Report progress after this many iterations. Set to Inf to turn it off
#' @return A list of treedater trees
#' @examples
#' \dontrun{
#' # make a random tree:
#' tre <- ape::rtree(50)
#' # sample times based on distance from root to tip:
#' sts <- setNames( ape::node.depth.edgelength( tre )[1:ape::Ntip(tre)], tre$tip.label)
#' # modify edge length to represent evolutionary distance with rate 1e-3:
#' tre$edge.length <- tre$edge.length * 1e-3
#' # treedater:
#' td <- dater( tre, sts =sts, clock='strict', s = 1000, omega0=.0015 )
#' gibbs_jitter( td )
#'}
#' @export
gibbs_jitter <- function(dtr, iter = 1e3, burn_pc = 20, returnTrees = 10, res = 100 , report = 10)
{
	#wishlist: version that modifies 1) tip dates 2) root height
	n <- length( dtr$tip )
	t <- c( dtr$sts[dtr$intree$tip.label] , dtr$Ti ) # current state
	sampleOrderNodes <- sample( (n+1):(n + dtr$Nnode),  replace=F) # order of nodes to sample

	td <- .make.tree.data (  dtr$intree, dtr$sts, dtr$s, cc = 10)
	td$minblen <- dtr$minblen #ugly

	nodes <- 1:(n + dtr$Nnode)
	node2edgei_list  <- lapply( nodes, function(x){
		which( dtr$intree$edge[,2] == x )
	})

	.sample.ti <- function(node )
	{
		# a sample/importance/resample algorithm with uniform proposal
		dgtrs <- td$daughters[node, ]
		a <- td$parent[ node ]
		if (any(is.na( c( dgtrs, a )))) return(NA)
		b1 <- td$tre$edge.length[ node2edgei_list[[dgtrs[1]]] ]
		b2 <- td$tre$edge.length[ node2edgei_list[[dgtrs[2]]] ]
		b3 <- td$tre$edge.length[ node2edgei_list[[ node  ]] ]

		tub <- min( t[dgtrs ] )
		tlb <- t[ a ]
		if ( tlb == tub ) return( NA )

		#tx <- seq( tlb, tub, l = 100 ) #TODO can probs do better than this
		tx <- runif( res , tlb , tub )

		# vectorised:
		u1s <-  t[ dgtrs[1] ] - tx
		u2s <-  t[ dgtrs[2] ] - tx
		u3s <-  tx-t[a]

		p1s <- dtr$theta * u1s / ( 1 + dtr$theta * u1s )
		p2s <- dtr$theta * u2s / ( 1 + dtr$theta * u2s )
		p3s <- dtr$theta * u3s / ( 1 + dtr$theta * u3s )

		if ( dtr$clock=='strict' ){
			lls <- dpois( round(b3 * dtr$s), u1s*dtr$mean.rate*dtr$s , log =T ) +
				dpois( round(b3 * dtr$s), u2s*dtr$mean.rate*dtr$s , log =T ) +
				dpois( round(b3 * dtr$s), u3s*dtr$mean.rate*dtr$s , log =T )
		} else if (dtr$clock == 'uncorrelated' ){
			lls <- dnbinom( round( b1 * dtr$s ) , dtr$r, 1 - p1s , log = TRUE ) +
				dnbinom( round(b2 * dtr$s), dtr$r , 1 - p2s, log =T ) +
				dnbinom( round(b3 * dtr$s), dtr$r, 1 - p3s , log =T )
		} else{
			stop('clock not implemented')
		}

		lls[is.na(lls)] <- -Inf
		if (max(lls)==-Inf) return(NA)
		w <- exp( lls - max( lls )  )
		if (sum(w)==0) {
			warning('All sample weights zero')
			return( NA )
		}
		tx [ sample(1:length(tx), size = 1, prob= w )]
	}

	X <- matrix( NA, nrow = length(t), ncol = iter)
	for (i in 1:iter){
		sampleOrderNodes <- sample( (n+1):(n + dtr$Nnode),  replace=F) # order of nodes to sample
		for ( node in sampleOrderNodes ){
			ti <- .sample.ti( node )
			if (!is.na( ti )) t[node] <- ti
		}
		X[, i] <- t

		if ( (i %% report)  == 0 ){
			print( paste( i, Sys.time() )  )
		}
	}

	# burn & sample t's
	ix <- round( seq( floor( burn_pc * iter/100), iter, l = returnTrees ) )
	X <- X[ , ix ]

	# return daters
	lapply( 1:ncol(X), function(i){
		t <- X[, i ]
		Ti <- t[ (n+1):(n + dtr$Nnode ) ]
		dtr$Ti <- Ti
		dtr$edge.length <- .Ti2blen(Ti, td )
		dtr
	})
}
