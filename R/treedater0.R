# gamma-poisson (NB) model for rate variation
#~ initial guess omega0 from rtt
#~ .omega0 -> Ti (assuming uniform omega_i)
#~ Ti -> optim r,p
#~ r,p -> optim omega_i
#~ omega_i -> Ti 
#~ repeat


require(ape)
require(mgcv)


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

.optim.r.gammatheta.nbinom0 <- function(  Ti, r0, gammatheta0, td)
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
		ov 
	}
	x0 <- c( lnr = unname(log(r0)), lngammatheta = unname( log(gammatheta0)))
	#o <- optim( par = x0, fn = of, method = 'BFGS' )
	o <- optim( par = x0, fn = of)
	r <- unname( exp( o$par['lnr'] ))
	gammatheta <- unname( exp( o$par['lngammatheta'] ))
	list( r = r, gammatheta=gammatheta, ll = -o$value)
}

.optim.Ti0 <- function( omegas, td , scale_var_by_rate){
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
		warning('Numerical error when performing least squares optimisation. Values may be approximate.')
		rv[is.na(rv)] <- max(rv, na.rm=T)
		rv <- .hack.times1(rv, td )
	}
	rv
}

# TODO could also try  package limSolve to add constraints T_i > 0	
.optim.Ti1 <- function( omegas, Ti, td ){
	blen <- .Ti2blen( Ti, td )
		A <- omegas * td$s *  td$A0 
		B <- td$B0 
		B[td$tipEdges] <- td$B0[td$tipEdges] -  unname( omegas[td$tipEdges] * td$sts2 )
		B <- B * td$s
		return( coef( lm ( B ~ A -1 , weights = 1 / (td$s*blen*omegas) ) ) ) # note blen in units of time
}
.Ti.ll1 <-  function( omegas, Ti, td ){
	blen <- .Ti2blen( Ti, td )
	A <- omegas * td$s *  td$A0 
	B <- td$B0 
	B[td$tipEdges] <- td$B0[td$tipEdges] -  unname( omegas[td$tipEdges] * td$sts2 )
	B <- B * td$s
	sum( dnorm( (B - A %*% Ti ), mean=0 , sd=sqrt((td$s*blen*omegas))  , log = T) )
}
# constraints using quad prog 
.optim.Ti2 <- function( omegas, td ){
		A <- omegas * td$A0 
		B <- td$B0
		B[td$tipEdges] <- td$B0[td$tipEdges] -  unname( omegas[td$tipEdges] * td$sts2 )
		
		# initial feasible parameter values:
		#p0 <- ( coef( lm ( B ~ A -1 , weights = td$W/omegas) ) )
		p0 <- ( coef( lm ( B ~ A -1 , weights = td$W) ) )
		if (any(is.na(p0))){
			warning('Numerical error when performing least squares optimisation. Values may be approximate.')
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
		o <- pcls(M)
	o
}


.optim.omegas.gammaPoisson1 <- function( Ti, r, gammatheta, td )
{
	blen <- .Ti2blen( Ti, td )
	o <- sapply( 1:nrow(td$tre$edge), function(k){
		sb <- (td$tre$edge.length[k]*td$s)
		lb <- qgamma(1e-5,  shape=r, scale = gammatheta*blen[k] )
		lam_star <- max(lb, gammatheta * blen[k] * (sb + r - 1) / (gammatheta * blen[k] + 1)  )
		ll <- 	dpois( max(0, round(sb)),lam_star, log=T )  + 
		 dgamma(lam_star, shape=r, scale = gammatheta*blen[k], log = T)
		c(lam_star / blen[k] / td$s, ll )
	}) 
	list( omegas = o[1,], ll = unname(sum( o[2,] )), lls  = unname(o[2,])  )
}

.optim.omega.poisson0 <- function(Ti, omega0, td)
{	
	blen <- .Ti2blen( Ti, td )
	
	of <- function(omega)
	{
		-sum( dpois( pmax(0, round(td$tre$edge.length*td$s)), td$s * blen * omega ,  log = T) )
	}
	o <- optimise(  of, lower = omega0 / 10, upper = omega0 * 10)
	list( omega = unname( o$minimum), ll = -unname(o$objective) )
}

.optim.omega.normal0 <- function(Ti, omega0, td)
{	
	blen <- .Ti2blen( Ti, td ) 
	blen <- pmax( quantile( blen, .05 ), blen ) 
	
	of <- function(omega)
	{
		-sum( dnorm( td$s * td$tre$edge.length,  sd = sqrt( td$s * blen * omega ) ,  log = T) )
	}
	o <- optimise(  of, lower = omega0 / 10, upper = omega0 * 10)
	list( omega = unname( o$minimum), ll = -unname(o$objective) )
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
		if (ub == tu & lb == tu) return(tu + td$minblen)
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

treedater = dater <- function(tre, sts, s=1e3
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
 , numStartConditions = 1
)
{ 
	# defaults
	CV_LB <- .07 # switch to poisson model below this value (coef of variation of gamma)
	scale_var_by_rate <- FALSE # better performance on lsd tests without this 
	cc <- 1
	
	numStartConditions <- max(0, round( numStartConditions )) # number of omega0 to try for optimisation
	
	EST_SAMP_TIMES <- TRUE
	EST_SAMP_TIMES_ERR <- 'estimateSampleTimes must specify a data frame with tip.label as row names and with columns `upper` and `lower`. You may also provide a named list of log density functions (improper priors for sample times).\n'
	if (is.null(estimateSampleTimes)) EST_SAMP_TIMES <- FALSE
	.estimateSampleTimes_densities <- estimateSampleTimes_densities
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
					estimateSampleTimes_densities[[tl]] <-  function(x,tl) dunif(x, min= estimateSampleTimes[tl,'lower'], max=estimateSampleTimes[tl,'upper'] , log = TRUE) #TODO maybe deprecate
				}
			}
		} else {
			stop(EST_SAMP_TIMES_ERR)
		}
		tiplabel_est_samp_times <- intersect( rownames(estimateSampleTimes), tre$tip.label)
		iedge_tiplabel_est_samp_times <- match( tiplabel_est_samp_times, tre$tip.label[tre$edge[,2]] )
	}
	
	if (is.null(names(sts))){
		names(sts) <- tre$tip.label
	}
	sts <- sts[tre$tip.label]
	
	intree_rooted <- TRUE
	if (!is.rooted(tre)){
		intree_rooted <- FALSE
		if (!quiet) cat( 'Tree is not rooted. Searching for best root position. Increase searchRoot to try harder.\n')
		searchRoot <- round( searchRoot )
		rtres <- .multi.rtt(tre, sts, topx=searchRoot)
		tds <- lapply( rtres, function(t) {
			#tryCatch( {
				dater( t, sts, s = s, omega0=omega0, minblen=minblen, maxit=maxit,abstol=abstol
				, strictClock = strictClock, temporalConstraints = temporalConstraints, quiet = quiet
				, estimateSampleTimes = estimateSampleTimes
				, estimateSampleTimes_densities = .estimateSampleTimes_densities  ) 
			#}, error = function(e) list( loglik = -Inf)) # 
		})
		lls <- sapply( tds, function(td) td$loglik )
		return ( tds [[ which.max( lls ) ]] )
	} else{
		if (!quiet) cat( 'Tree is rooted. Not estimating root position.\n')
	}
	if (is.na(minblen)){
		minblen <- diff(range(sts))/100/ length(sts) #TODO choice of this parm is difficult, may require sep optim / crossval
	}
	if (!is.na(omega0) & numStartConditions > 0 ){
		warning('omega0 provided incompatible with numStartConditions > 0. Setting numStartConditions to zero.')
		numStartConditions <- 0
	}
	if (is.na(omega0)){
		# guess
		#omega0 <- estimate.mu( tre, sts )
		g0 <- lm(ape::node.depth.edgelength(tre)[1:length(sts)] ~ sts, na.action = na.omit)
		omega0sd <- summary( g0 )$coef[2,2]
		omega0 <- unname( coef(g0)[2] )
		omega0s <- qnorm( unique(sort(c(.5, seq(.025, .975, l=numStartConditions*2) )))  , omega0, sd = omega0sd )
		omega0s <- omega0s[ omega0s > 0 ]
		if (!quiet){
			cat('initial rates:\n')
			print(omega0s)
		}
	} else{
		omega0s <- c( omega0 )
	}
	td <- .make.tree.data(tre, sts, s, cc )
	td$minblen <- minblen
	
	#lapply( omega0s , function(omega0){
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
				Ti <- .optim.Ti2( omegas, td)
			} else{
				Ti <- .optim.Ti0( omegas, td, scale_var_by_rate )
			}
			
			if ( (1 / sqrt(r)) < CV_LB){
				# switch to poisson model
				o <- .optim.omega.poisson0(Ti, .mean.rate(Ti, r, gammatheta, omegas, td), td)
				#o <- .optim.omega.normal0(Ti, .mean.rate(Ti, r, gammatheta, omegas, td), td)
				gammatheta <- unname(o$omega)
				if (!is.infinite(r)) lastll <- -Inf # the first time it switches, do not do likelihood comparison 
				r <- Inf#unname(o$omega)
				ll <- o$ll
				edge_lls <- 0
				omegas <- rep( gammatheta, length(omegas))
			} else{
				o <- .optim.r.gammatheta.nbinom0(  Ti, r, gammatheta, td)
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
				cat('iter, omegas, T, r, theta, logLik\n')
				print(c( iter))
				print(summary(omegas))
				print(summary(Ti))
				print( r)
				print( gammatheta)
				print( ll)
			}
			if ( ll >= lastll ){
				rv <- list( omegas = omegas, r = unname(r), theta = unname(gammatheta), Ti = Ti
				 , meanRate = .mean.rate(Ti, r, gammatheta, omegas, td)
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
			omega2ll <- ll
		}
		#rv
	}
	#) -> rvs 
	
	#~ 	rvs_lls <- sapply( rvs, function(rv) rv$loglik )
	#~ 	rv <- rvs[[which.max( rvs_lls )]]
	
	rv <- bestrv
	
	.tre <- tre
	td$minblen <- -Inf; blen <- .Ti2blen( rv$Ti, td )
	tre$edge.length <- blen 
	
	rv$tre <- tre
	rv$edge <- tre$edge
	rv$edge.length <- tre$edge.length
	rv$tip.label <- tre$tip.label
	rv$Nnode <- tre$Nnode
	
	rv$timeOfMRCA <- min(rv$Ti)
	rv$timeToMRCA <- max(sts) - rv$timeOfMRCA
	rv$s <- s
	rv$sts <- sts
	rv$minblen <- minblen
	rv$intree <- .tre 
	rv$coef_of_variation <- 1 / sqrt(rv$r)
	rv$clock <- ifelse( is.infinite(rv$r), 'strict', 'relaxed')
	rv$intree_rooted <- intree_rooted
	rv$temporalConstraints <- temporalConstraints
	rv$estimateSampleTimes <- estimateSampleTimes
	rv$EST_SAMP_TIMES = EST_SAMP_TIMES
	if (!EST_SAMP_TIMES) rv$estimateSampleTimes <- NULL
	rv$estimateSampleTimes_densities <- estimateSampleTimes_densities
	class(rv) <- c('treedater', 'phylo')
	rv
}

print.treedater <- function(x, ...){
    cl <- oldClass(x)
    oldClass(x) <- cl[cl != "treedater"]
    print(x$tre)
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

summary.treedater <- function(x) {
    stopifnot(inherits(x, "treedater"))
    print.treedater( x )
}
