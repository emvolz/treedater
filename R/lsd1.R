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
	
	W <- abs( 1/((tre$edge.length + cc / s)/s) ) 
	list( A0 = A, B0 = B, W0 = W, n = n,  tipEdges=tipEdges
	 , i_tip_edge2label = i_tip_edge2label
	 , sts2 = sts 
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
		ps <- pmin(.999, gammatheta*blen / ( 1+ gammatheta * blen ) )
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
		if (scale_var_by_rate){
			return( coef( lm ( B ~ A -1 , weights = td$W/omegas) ) )
		} else{
			return( coef( lm ( B ~ A -1 , weights = td$W) ) ) 
		}
}

# TODO could also try  package limSolve to add constraints T_i > 0	
.optim.Ti1 <- function( omegas, Ti, td ){
	blen <- .Ti2blen( Ti, td )
		A <- omegas * td$s *  td$A0 
		B <- td$B0 
		B[td$tipEdges] <- td$B0[td$tipEdges] -  unname( omegas[td$tipEdges] * td$sts2 )
		B <- B * td$s
		return( coef( lm ( B ~ A -1 , weights = 1 / (td$s*blen*omegas) ) ) )
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
		p0 <- ( coef( lm ( B ~ A -1 , weights = td$W/omegas) ) )
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
			,w=td$W/omegas
		)
		o <- pcls(M)
	o
}

.optim.omegas.gammaPoisson0 <- function( Ti, r, gammatheta, td )
{
	blen <- .Ti2blen( Ti, td )
	
	omegas <- rep(NA, length(blen))
	
	for (k in 1:nrow(td$tre$edge) ){
		#lb <- 1e-3  
		lb <- qgamma(1e-5,  shape=r, scale = gammatheta*blen[k] )
		ub <- qgamma(1-1e-5,  shape=r, scale = gammatheta*blen[k] )
		of <- function(lam_k){
			-dpois( max(0, round(td$tre$edge.length[k]*td$s)),lam_k, log=T )  - 
			 dgamma(lam_k, shape=r, scale = gammatheta*blen[k], log = T)
		}
		lam_k <- optimise( of, lower = lb, upper = ub)$minimum
		omegas[k] <- lam_k / blen[k] / td$s
	}
	omegas
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


.mean.rate <- function(Ti, r, gammatheta, omegas, td)
{
	if (gammatheta<=0){
		# this is a poisson model
		return(r)
	}	
	blen <- .Ti2blen( Ti, td )
	
	sum( omegas * blen ) / sum(blen)
}

.hack.times1 <- function(Ti, td)
{
	t <- c( td$sts2[td$tre$tip.label], Ti)
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
 , abstol = .01
 , searchRoot = 5
 , quiet = FALSE
 , temporalConstraints = TRUE
 , strictClock = FALSE
)
{ 
	THETA_LB <- 1e-3
	scale_var_by_rate <- TRUE
	cc <- 10
	if (is.null(names(sts))){
		names(sts) <- tre$tip.label
	}
	sts <- sts[tre$tip.label]
	if (!is.rooted(tre)){
		if (!quiet) cat( 'Tree is not rooted. Searching for best root position. Increase searchRoot to try harder.\n')
		searchRoot <- round(searchRoot )
		rtres <- .multi.rtt(tre, sts, topx=searchRoot)
		tds <- lapply( rtres, function(t) dater( t, sts, s = s, omega0=omega0, minblen=minblen, maxit=maxit,abstol=abstol
		 , strictClock = strictClock, temporalConstraints = temporalConstraints, quiet = quiet) )
		lls <- sapply( tds, function(td) td$loglik )
		return ( tds [[ which.max( lls ) ]] )
	} else{
		if (!quiet) cat( 'Tree is rooted. Not estimating root position.\n')
	}
	sts <- sts[tre$tip.label]
	if (is.na(minblen)){
		minblen <- (max(sts) - min(sts)) / length(sts) #TODO choice of this parm is difficult, may require sep optim / crossval
	}
	if (is.na(omega0)){
		# guess
		n <- length( tre$tip.label)
		d2root <- setNames(  dist.nodes( tre )[(n+1),1:n], tre$tip.label)
		omega0 <- coef( lm( d2root ~ sts) )[2]
		if (!quiet){
			cat('initial rate:\n')
			print(omega0)
		}
	}
	td <- .make.tree.data(tre, sts, s, cc )
	td$minblen <- minblen
	
	# initial gamma parms with small variance
	r0 = r <- ifelse( strictClock, omega0 , omega0 * td$s / 1e-2)
	gammatheta0 = gammatheta <- ifelse( strictClock, 0,  1e-2 )
	
	done <- FALSE
	lastll <- -Inf
	iter <- 0
	nEdges <- nrow(tre$edge)
	omegas <- rep( omega0,  nEdges )
	.trace <- list()
	while(!done){
		if (temporalConstraints){
			Ti <- .optim.Ti2( omegas, td)
		} else{
			Ti <- .optim.Ti0( omegas, td, scale_var_by_rate )
		}
		
		if (gammatheta < THETA_LB){
			# switch to poisson model
			# r is now rate parm
			o <- .optim.omega.poisson0(Ti, .mean.rate(Ti, r, gammatheta, omegas, td), td)
			gammatheta <- 0
			r <- unname(o$omega)
			ll <- o$ll
			omegas <- rep( r, length(omegas))
		} else{
			o <- .optim.r.gammatheta.nbinom0(  Ti, r, gammatheta, td)
			r <- o$r
			ll <- o$ll
			gammatheta <- o$gammatheta
			omegas <- .optim.omegas.gammaPoisson0( Ti, o$r, o$gammatheta, td ) 
		}
		
		ti_ll <- .Ti.ll1( omegas, Ti, td ) #
		ll <- ll + ti_ll
		
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
		.trace[[length(.trace)+1]] <- list( omegas = omegas, r = r, theta = gammatheta, Ti = Ti
		 , meanRate = .mean.rate(Ti, r, gammatheta, omegas, td)
		 , loglik = ll )
		
		# check convergence
		iter <- iter + 1
		if (iter > maxit) done <- TRUE
		
		if ( abs( ll - lastll ) < abstol) done <- TRUE
		if (ll < lastll) done <- TRUE
		
		lastll <- ll
	}
	i <- which.max( sapply( .trace, function(x) x$loglik) )
	rv <- .trace[[i]]
	td$minblen <- -Inf; blen <- .Ti2blen( rv$Ti, td )
	tre$edge.length <- blen 
	#rv$trace <- .trace
	rv$tre <- tre
	rv$timeOfMRCA <- min(rv$Ti)
	rv$timeToMRCA <- max(sts) - rv$timeOfMRCA
	rv
}
