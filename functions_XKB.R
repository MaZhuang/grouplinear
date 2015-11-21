# XKB (2012) code

library(isotone)

# bayes rule for fixed lambda,mu
thetahat <- function(X,A,lambda,mu){
	lambda/(lambda+A) * X + A/(lambda+A) * mu
}

# d/dlambda{SURE(lambda,mu=mu.hat.SURE(lambda))} (proportional to)
g <- function(lambda,X,A){  
	sum( A^2/(lambda+A)^3 * (X-(  sum( A^2/(lambda+A)^2 * X ) / sum( A^2/(lambda+A)^2 )  ))^2 - A^2/(lambda+A)^2 )
	#equivalent to the following(which is just easier to read):
	#mu <- sum( A^2/(lambda+A)^2 * X ) / sum( A^2/(lambda+A)^2 )
	#sum( A^2/(lambda+A)^3 * (X-mu)^2 - A^2/(lambda+A)^2 )
}

# SURE(lambda,mu)
f <- function(par,X,A){  
	lambda <- par[1]
	mu <- par[2]
	sum(  A/(lambda+A)^2 * ( A * (X-mu)^2 + lambda^2 - A^2 )  )
}

# SURE.G(lambda)
f.G <- function(lambda,X,A){  
	sum(  ( A/(A+lambda) )^2 * (X - mean(X))^2 + A/(A+lambda) * (lambda - A + 2/p * A)  )
}


thetahat.M <- function(X,A){
	lambda.sure <- ifelse( g(0,X=X,A=A)*g(max(A)*1000,X=X,A=A) < 0, uniroot(g,c(0,max(A)*1000),X=X,A=A, tol=1e-9)$root, optim(c(mean(   pmax(  ( X-mean(X) )^2 - 	A,0  )   ), mean(X)),f,X=X,A=A, method = "L-BFGS-B",lower=c(0,-Inf))$par[1] )
	mu.sure <- sum( A^2/(lambda.sure+A)^2 * X ) / sum( A^2/(lambda.sure+A)^2 )
	thetahat(X,A,lambda.sure,mu.sure)
}

thetahat.G <- function(X,A){
lambda <- optimize(f.G,lower=0,upper=1000,X=X,A=A)$minimum
thetahat(X,A,lambda,mean(X))
}

thetahat.SG <- function(X,A){
	p <- length(X)
	fit <- gpava( z = A, y = A * (1-1/p) / (X-mean(X))^2, weights = (X-mean(X))^2, solver = weighted.mean, ties="primary" )
	bhat <- pmin(  pmax( fit$x,0 ),1  )
	(1-bhat) * X + bhat * mean(X)
}

#example
p <- 100
A <- runif(p,.1,1)
theta <- rnorm(p)
X <- rnorm(p,theta,sqrt(A))
thetahat.M(X,A)
thetahat.G(X,A)
thetahat.SG(X,A)


# plot(X,thetahat.M(X,A), cex=.5, pch=16, ylim = range(thetahat.M(X,A), thetahat.SG(X,A)))
# points(X,thetahat.SG(X,A), cex=.5, pch=16, col='blue')


# SURE^G(lambda)

sure.G <- function(lambda,x,v){
  x.bar <- mean(x)
  mean(  v^2/(v+lambda)^2 * (x-x.bar)^2 + v/(v+lambda) * (lambda - v + 2/n * v)  )  
}

# lambda.vec <- ppoints(100)*100
# sure.G(lambda.vec,X,A)
# y <- sapply(lambda.vec,sure.G,x=x,v=v)
# plot(lambda.vec,y,pch=16,cex=.5)




