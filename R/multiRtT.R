###
# based on ape::rtt developed by  Rosemary McCloskey, BC Centre for Excellence in HIV/AIDS, released under GPL
# retrieved December 2016 from https://github.com/cran/ape/blob/master/R/rtt.R
# returns multiple rooted trees reflecting best fits of root-to-tip regression
.multi.rtt <-
function (t, tip.dates, topx=1, ncpu = 1, objective = "correlation",  opt.tol = .Machine$double.eps^0.25) 
{
	topx <- max( 1, topx)
	if (objective == "correlation") 
		objective <- function(x, y) cor.test(y, x)$estimate
	else if (objective == "rsquared") 
		objective <- function(x, y) summary(lm(y ~ x))$r.squared
	else if (objective == "rms") 
		objective <- function(x, y) -summary(lm(y ~ x))$sigma^2
	else stop("objective must be one of \"correlation\", \"rsquared\", or \"rms\"")
	ut <- unroot(t)
	dist <- dist.nodes(ut)[, 1:(ut$Nnode + 2)]
	f <- function(x, parent, child) {
		edge.dist <- x * dist[parent, ] + (1 - x) * dist[child, 
			]
		objective(tip.dates, edge.dist)
	}
	obj.edge <- if (ncpu > 1) 
		unlist(parallel::mclapply(1:nrow(ut$edge), function(e) {
			opt.fun <- function(x) f(x, ut$edge[e, 1], ut$edge[e, 
				2])
			optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$objective
		}, mc.cores = ncpu))
	else apply(ut$edge, 1, function(e) {
		opt.fun <- function(x) f(x, e[1], e[2])
		optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$objective
	})
	obj.edge[ ut$edge.length<=0 ]  <- -Inf # excludes nodes with zero branch parents 
	
	best.edges <- order( obj.edge, decreasing=TRUE)[1:topx]
	#for (best.edge in best.edges )
	lapply(best.edges , function(best.edge){
		best.edge.parent <- ut$edge[best.edge, 1]
		best.edge.child <- ut$edge[best.edge, 2]
		best.edge.length <- ut$edge.length[best.edge]
		opt.fun <- function(x) f(x, best.edge.parent, best.edge.child)
		best.pos <- optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$maximum
		new.root <- list(edge = matrix(c(2L, 1L), 1, 2), tip.label = "new.root", 
			edge.length = 1, Nnode = 1L, root.edge = 1)
		class(new.root) <- "phylo"
		ut2 <- bind.tree(ut, new.root, where = best.edge.child, position = best.pos * 
			best.edge.length)
		ut2 <- collapse.singles(ut2)
		ut2 <- root(ut2, "new.root")
		x <- drop.tip(ut2, "new.root")
		if (!is.rooted(x)) return(NULL)
		x
	})-> tres
	tres[ !sapply( tres, is.null) ] 
}
