# Group-linear Functions

# "2014-11-13 11:33:34 EST"

## spherically symmetric estimator with c_n = c^*_n
spher <- function(x.,v.){
n. <- length(x.)
if ( (n.==1) | (var(x.)==0) ) x. else {
	cstar <- max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0)
	bhat <- min( cstar*mean(v.)/var(x.), 1 )
	x. - bhat*(x. - mean(x.))
	}
}


## spherically symmetric estimator with c_n = c^*_n, shrinkage toward zero
spher.zero <- function(x.,v.){
  n. <- length(x.)
  cstar <- max( 1-2*( max(v.)/mean(v.) )/n., 0)
  bhat <- min( cstar*mean(v.)/mean(x.^2), 1 )
  (1- bhat)*x.
}

## function that returns the common bhat (replicated)
spher.bhat <- function(x.,v.){
  n. <- length(x.)
  if ( (n.==1) | (var(x.)==0) ) x. else {
    cstar <- max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0)
    bhat <- min( cstar*mean(v.)/var(x.), 1 )
    return(rep(bhat,n.))
  }
}

## group-linear estimator

grouplinear <- function( x,v,nbreak=floor(length(x)^(1/3)) ){  # default: bin log(v) into same NUMBER (=n^(1/3) of intervals
  n <- length(x)
  splitby=cut(log(v),breaks=nbreak, labels=F)
	xsub <- split(x,splitby)
	vsub <- split(v,splitby)
	indexsub <- split(1:n,splitby)
	thetahatsub <- mapply(spher,xsub,vsub)
	indexsub.unlist <- as.vector( unlist(indexsub) )
	thetahatsub.unlist <- as.vector( unlist(thetahatsub) )
	thetahat <- thetahatsub.unlist[order(indexsub.unlist)]	
	return(thetahat)
}

## group-linear estimator with shrinkage toward zero

grouplinear.zero <- function( x,v,nbreak=floor(length(x)^(1/3)) ){  # default: bin log(v) into same NUMBER (=n^(1/3) of intervals
  n <- length(x)
  splitby=cut(log(v),breaks=nbreak, labels=F)
  xsub <- split(x,splitby)
  vsub <- split(v,splitby)
  indexsub <- split(1:n,splitby)
  thetahatsub <- mapply(spher.zero,xsub,vsub)
  indexsub.unlist <- as.vector( unlist(indexsub) )
  thetahatsub.unlist <- as.vector( unlist(thetahatsub) )
  thetahat <- thetahatsub.unlist[order(indexsub.unlist)]	
  return(thetahat)
}

# # ## example
# n <- 300
# v <- runif(n,.1,1)
# theta <- v-mean(v)
# x <- rnorm(n,theta,sd=sqrt(v))
# grouplinear(x,v)
# mean( (grouplinear(x,v)-theta)^2 )   
# mean( (grouplinear.zero(x,v)-theta)^2 )   


## sure for grouplinear estimator

sure.spher <- function(x.,v.){
  n. <- length(x.)
 # cstar <- max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0) ##modified
  if (n.==0) {NULL 
  } else if ( (n.<3) ) {sum(v.)  #| (var(x.)==0) 
  }
  else if (max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0)==0){
  	sum(v.) 
  }
  else if (var(x.)==0){
  	(2-n.)/n.*sum(v.)+sum((x.-mean(x.))^2)
  }
  else {	# can set sure to an arbitrary value if var(x.)=0, since this event is of measure zero
	cstar <- max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0) ##modified
	b <- cstar * mean(v.)/var(x.)
    b <- min(1,b)
	db <- -cstar * mean(v.)/(var(x.))^2 * as.numeric( cstar * mean(v.)/var(x.) < 1 )##
	sum(   v. + ( b * (x.-mean(x.)) )^2 - 2 * v. * (  (1-1/n.) * b + 2 * (x.-mean(x.))^2 * db/(n.-1)  )   )
	}
}



sure.spher.zero <- function(x.,v.){
  if (n.==0) {0 
  }else if ( (n.==1) | (var(x.)==0) ) {sum(v.) 
  }else {	# can set sure to an arbitrary value if var(x.)=0, since this event is of measure zero
    cstar <- 1-2*max(v.)/sum(v.) * ( 1-2*max(v.)/sum(v.) > 0 )
    b <- cstar * mean(v.)/mean(x.^2)
    b <- min(1,b)
    db <- -cstar * mean(v.)/mean(x.^2)^2 * ( cstar * mean(v.)/mean(x.^2) < 1 )
    sum(   v. + ( b * (x.) )^2 - 2 * v. * (  b + 2 * x.^2 * db/n.  )   )
  }
}

sure.grouplinear <- function(x,v,nbreak){ #nbreak=num of bins
	n <- length(x)
	splitby=cut(log(v),breaks=nbreak)
	xsub <- split(x,splitby)
	vsub <- split(v,splitby)
	suresub <- mapply(sure.spher,xsub,vsub)   #modified
	sum(suresub)/n
}

sure.grouplinear.zero <- function(x,v,nbreak){ #nbreak=num of bins
  n <- length(x)
  splitby=cut(log(v),breaks=nbreak)
  xsub <- split(x,splitby)
  vsub <- split(v,splitby)
  suresub <- mapply(sure.spher.zero,xsub,vsub)
  sum(suresub)/n
}


