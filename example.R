# Baseball

if(!exists("foo", mode="function")) source("functions.R")
if(!exists("foo", mode="function")) source("functions_XKB.R")
if(!exists("foo", mode="function")) source("dynamic_sure.R")
path=getwd()
datapath=paste(path,'Brown_batting_data.txt',sep='/')
bat.raw <- read.table(datapath, header=TRUE, sep=",", quote="")

# prepare the data for analysis
#bat.raw <- read.table("~/desktop/Example/Brown_batting_data.txt", header=TRUE, sep=",", quote="")

bat <- bat.raw
bat$N1 <- bat$AB.4. + bat$AB.5. + bat$AB.6.  # total number at-bats for 1st period
bat$N2 <- bat$AB.7. + bat$AB.8. + bat$AB.9.10.  # total number at-bats for 2nd period
bat$H1 <- bat$H.4. + bat$H.5. + bat$H.6.  # total number hits for 1st period
bat$H2 <- bat$H.7. + bat$H.8. + bat$H.9.10.  # total number hits for 2nd period
bat$R1 <- bat$H1/bat$N1  # batting avg for 1st period
bat$R2 <- bat$H2/bat$N2  # batting avg for 2nd period
bat$X1 <- asin(  sqrt( (bat$H1+1/4)/(bat$N1+1/2) )  )  # transformed batting avg for 1st period
bat$X2 <- asin(  sqrt( (bat$H2+1/4)/(bat$N2+1/2) )  )  # transformed batting avg for 2nd period
bat <- bat[bat$N1 > 10, c('First.Name','Last.Name','Pitcher.','N1','N2','H1','H2','X1','X2')]  # keep only records with N1>=11


index=order(bat$N1,decreasing=TRUE)
bat=bat[index,]
#bat=bat[1:n,]

ind <- bat$N2>10  # indicator for records with N2>=11 (among those with N1>=11)

tse.hat.zero <- sum(   (  ( bat$X2 - bat$X1 )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )

# grand mean
tse.hat.delta.gm <- sum(   (  ( bat$X2 - mean(bat$X1) )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.hat.delta.gm/tse.hat.zero

# XKB theta.hat.M
delta.M <- thetahat.M(bat$X1,1/(4 * bat$N1))
tse.hat.delta.M <- sum(   (  ( bat$X2 - delta.M )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.hat.delta.M/tse.hat.zero

# XKB theta.hat.SG
delta.SG <- thetahat.SG(X=bat$X1,A=1/(4 * bat$N1))
tse.hat.delta.SG <- sum(   (  ( bat$X2 - delta.SG )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.hat.delta.SG/tse.hat.zero

# group-linear
  # num bins = n^1/3
delta.gl <- grouplinear(x=bat$X1, v=1/(4 * bat$N1))
tse.hat.delta.gl <- sum(   (  ( bat$X2 - delta.gl )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.hat.delta.gl/tse.hat.zero



#dynamic
delta.dynamic=GroupSure(bat$X1,1/(4 * bat$N1))
tse.hat.delta.dynamic <- sum(   (  ( bat$X2 - delta.dynamic )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.hat.delta.dynamic/tse.hat.zero




  # URE
  # i) split into k intervals of equal length on log(v)
min.diff <- min(diff( sort(log( 1/(4 * bat$N1) )) )[diff( sort(log( 1/(4 * bat$N1) )) )>0])  # min_{i,j: v_i != v_j} |v_i-v_j|
kmax <- ceiling( diff(range(log( 1/(4 * bat$N1) )))/min.diff )
sure.vec <- rep(NA,kmax)

sure.vec[1] <- sure.spher(bat$X1,1/(4 * bat$N1))
for (k in 2:30){
	sure.vec[k] <- sure.grouplinear(bat$X1,1/(4 * bat$N1),nbreak=k)
}

khat.sure <- which.min(sure.vec) # 4
delta.gl.sure <- if(khat.sure>1) grouplinear( bat$X1,1/(4 * bat$N1),nbreak=khat.sure ) else spher(bat$X1,1/(4 * bat$N1))
tse.hat.delta.gl.sure <- sum(   (  ( bat$X2 - delta.gl.sure )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.hat.delta.gl.sure/tse.hat.zero  # 0.34











  # oracle
rel.tse.breaks <- rep(NA,20)
delta.gl <- spher(x=bat$X1, v=1/(4 * bat$N1))
tse.hat.delta.gl <- sum(   (  ( bat$X2 - delta.gl )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
rel.tse.breaks[1] <- tse.hat.delta.gl/tse.hat.zero
for(i in 2:20){
  delta.gl <- grouplinear(x=bat$X1, v=1/(4 * bat$N1),nbreak = i)
  tse.hat.delta.gl <- sum(   (  ( bat$X2 - delta.gl )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  rel.tse.breaks[i] <- tse.hat.delta.gl/tse.hat.zero
}

tse.hat.delta.gl.ol <- min(rel.tse.breaks) #.284
k.ol <- which.min(rel.tse.breaks) # 10








# Analysis with covariates (Using subspaces shrinkage as in Tan's ms)
  # 1 + AB + pitcher
# head( cbind(bat$Pitcher.,bat$N1) )
lm1 <- lm(bat$X1 ~ as.factor(bat$Pitcher.) + bat$N1, weights = 4 * bat$N1)
y <- resid(lm1)
sig.hat <- summary(lm1)$sigma
  #compute covariance of the residual vector
cov.y <- diag(1/(4 * bat$N1)) - model.matrix(lm1) %*% (vcov(lm1)/sig.hat^2) %*% t(model.matrix(lm1)) # Note: division by sig.hat^2 makes a difference!
s <- svd(cov.y)
#sort(s$d)[1:5]
u1 <- (s$u)[,1:564]
  # transform (to a vector with full-rank diagonal cov matrix)
y.white <- crossprod(u1,y)
#estimate.white <- grouplinear(x=y.white, v=s$d[1:564])
estimate.white <- grouplinear.zero(x=y.white, v=s$d[1:564])
  # back-transform
estimate <- u1 %*% estimate.white
  # final estimate
delta.gl.sub <- fitted(lm1) + estimate

tse.hat.delta.gl.sub <- sum(   (  ( bat$X2 - delta.gl.sub )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.hat.delta.gl.sub/tse.hat.zero

# 1 + AB * pitcher
# head( cbind(bat$Pitcher.,bat$N1) )
lm1 <- lm(bat$X1 ~ as.factor(bat$Pitcher.) * bat$N1, weights = 4 * bat$N1)
y <- resid(lm1)
sig.hat <- summary(lm1)$sigma
#compute covariance of the residual vector
cov.y <- diag(1/(4 * bat$N1)) - model.matrix(lm1) %*% (vcov(lm1)/sig.hat^2) %*% t(model.matrix(lm1)) # Note: division by sig.hat^2 makes a difference!
s <- svd(cov.y)
#sort(s$d)[1:5]
u1 <- (s$u)[,1:563]
# transform (to a vector with full-rank diagonal cov matrix)
y.white <- crossprod(u1,y)
#estimate.white <- grouplinear(x=y.white, v=s$d[1:563])
estimate.white <- grouplinear.zero(x=y.white, v=s$d[1:563])
# back-transform
estimate <- u1 %*% estimate.white
# final estimate
delta.gl.sub <- fitted(lm1) + estimate

tse.hat.delta.gl.sub <- sum(   (  ( bat$X2 - delta.gl.sub )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.hat.delta.gl.sub/tse.hat.zero

# num bins = 2:10
rel.tse.breaks <- rep(NA,20)
for(i in 2:20){
  estimate.white <- grouplinear.zero(x=y.white, v=s$d[1:563],nbreak = i)
  # back-transform
  estimate <- u1 %*% estimate.white
  # final estimate
  delta.gl.sub <- fitted(lm1) + estimate
  tse.hat.delta.gl.sub <- sum(   (  ( bat$X2 - delta.gl.sub )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  rel.tse.breaks[i] <- tse.hat.delta.gl.sub/tse.hat.zero
}
# compute rel.tse for a single bin
  estimate.white <- spher.zero(x.=y.white, v.=s$d[1:563])
# back-transform
estimate <- u1 %*% estimate.white
# final estimate
delta.gl.sub <- fitted(lm1) + estimate
tse.hat.delta.gl.sub <- sum(   (  ( bat$X2 - delta.gl.sub )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
rel.tse.breaks[1] <- tse.hat.delta.gl.sub/tse.hat.zero

cbind(nbreaks=order(rel.tse.breaks),rel.tse=sort(rel.tse.breaks))

# nbreaks   rel.tse
# [1,]       1 0.1750941
# [2,]      16 0.1793561
# [3,]       3 0.1841942
# [4,]       7 0.1875441
# [5,]       5 0.1882745
# [6,]      13 0.1920799
# [7,]      17 0.1934835
# [8,]      14 0.1942529
# [9,]      18 0.1944820
# [10,]      10 0.1949698
# [11,]       9 0.1955755
# [12,]       2 0.1969034
# [13,]      15 0.1972869
# [14,]      19 0.1976090
# [15,]       4 0.2006492
# [16,]       8 0.2009908
# [17,]      11 0.2068433
# [18,]       6 0.2075451
# [19,]      12 0.2082239
# [20,]      20 0.2083171

# SURE-binning (1+AB*pitcher)
sure.vec <- rep(NA,30)
sure.vec[1] <- sure.spher.zero(x=y.white, v=s$d[1:564])
for (k in 2:30){
  sure.vec[k] <- sure.grouplinear.zero(x=y.white, v=s$d[1:564],nbreak=k)
}
( khat.sure <- which.min(sure.vec) )
# delta.gl.sure <- grouplinear( bat$X1,1/(4 * bat$N1),nbreak=khat.sure)
# tse.hat.delta.gl.sure <- sum(   (  ( bat$X2 - delta.gl.sure )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
# tse.hat.delta.gl.sure/tse.hat.zero  # 0.34



  # 1 + pitcher
# head( cbind(bat$Pitcher.,bat$N1) )
lm1 <- lm(bat$X1 ~ as.factor(bat$Pitcher.), weights = 4 * bat$N1)
y <- resid(lm1)
sig.hat <- summary(lm1)$sigma
#compute covariance of the residual vector
cov.y <- diag(1/(4 * bat$N1)) - model.matrix(lm1) %*% (vcov(lm1)/sig.hat^2) %*% t(model.matrix(lm1)) # Note: division by sig.hat^2 makes a difference!
s <- svd(cov.y)
#sort(s$d)[1:5]
u1 <- (s$u)[,1:565]
# transform (to a vector with full-rank diagonal cov matrix)
y.white <- crossprod(u1,y)
#estimate.white <- grouplinear(x=y.white, v=s$d[1:565])
estimate.white <- grouplinear.zero(x=y.white, v=s$d[1:565])
# back-transform
estimate <- u1 %*% estimate.white
# final estimate
delta.gl.sub <- fitted(lm1) + estimate

tse.hat.delta.gl.sub <- sum(   (  ( bat$X2 - delta.gl.sub )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.hat.delta.gl.sub/tse.hat.zero

# num bins = 2:10
rel.tse.breaks <- rep(NA,20)
for(i in 2:20){
  estimate.white <- grouplinear.zero(x=y.white, v=s$d[1:565],nbreak = i)
  # back-transform
  estimate <- u1 %*% estimate.white
  # final estimate
  delta.gl.sub <- fitted(lm1) + estimate
  tse.hat.delta.gl.sub <- sum(   (  ( bat$X2 - delta.gl.sub )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  rel.tse.breaks[i] <- tse.hat.delta.gl.sub/tse.hat.zero
}
# compute rel.tse for a single bin
estimate.white <- spher.zero(x.=y.white, v.=s$d[1:565])
# back-transform
estimate <- u1 %*% estimate.white
# final estimate
delta.gl.sub <- fitted(lm1) + estimate
tse.hat.delta.gl.sub <- sum(   (  ( bat$X2 - delta.gl.sub )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
rel.tse.breaks[1] <- tse.hat.delta.gl.sub/tse.hat.zero

cbind(nbreaks=order(rel.tse.breaks),rel.tse=sort(rel.tse.breaks))




## pitchers only

bat.pitch <- bat[bat$Pitcher.==1,]

# estimating TSE for various estimators
# run: functions.R(current folder), functions_XKB.R

ind.pitch <- bat.pitch$N2>10  # indicator for records with N2>=11 (among those with N1>=11)

tse.hat.zero <- sum(   (  ( bat.pitch$X2 - bat.pitch$X1 )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )

# grand mean
tse.hat.delta.gm <- sum(   (  ( bat.pitch$X2 - mean(bat.pitch$X1) )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
tse.hat.delta.gm/tse.hat.zero

# XKB theta.hat.M
delta.M <- thetahat.M(bat.pitch$X1,1/(4 * bat.pitch$N1))
tse.hat.delta.M <- sum(   (  ( bat.pitch$X2 - delta.M )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
tse.hat.delta.M/tse.hat.zero

# XKB theta.hat.SG
delta.SG <- thetahat.SG(bat.pitch$X1,1/(4 * bat.pitch$N1))
tse.hat.delta.SG <- sum(   (  ( bat.pitch$X2 - delta.SG )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
tse.hat.delta.SG/tse.hat.zero

# group-linear 
  # n^-1/3 bins
delta.gl <- grouplinear(x=bat.pitch$X1, v=1/(4 * bat.pitch$N1))
tse.hat.delta.gl <- sum(   (  ( bat.pitch$X2 - delta.gl )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
tse.hat.delta.gl/tse.hat.zero
  # oracle
rel.tse.breaks <- rep(NA,20)
delta.gl <- spher(x=bat.pitch$X1, v=1/(4 * bat.pitch$N1))
tse.hat.delta.gl <- sum(   (  ( bat.pitch$X2 - delta.gl )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
rel.tse.breaks[1] <- tse.hat.delta.gl/tse.hat.zero
for(i in 2:20){
  delta.gl <- grouplinear(x=bat.pitch$X1, v=1/(4 * bat.pitch$N1),nbreak = i)
  tse.hat.delta.gl <- sum(   (  ( bat.pitch$X2 - delta.gl )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
  rel.tse.breaks[i] <- tse.hat.delta.gl/tse.hat.zero
}

tse.hat.delta.gl.ol <- min(rel.tse.breaks) #.168
k.ol <- which.min(rel.tse.breaks) #1


# URE
# i) split into k intervals of equal length on log(v)
min.diff <- min(diff( sort(log( 1/(4 * bat.pitch$N1) )) )[diff( sort(log( 1/(4 * bat.pitch$N1) )) )>0])  # min_{i,j: v_i != v_j} |v_i-v_j|
kmax <- ceiling( diff(range(log( 1/(4 * bat.pitch$N1) )))/min.diff )
sure.vec <- rep(NA,kmax)

sure.vec[1] <- sure.spher(bat.pitch$X1,1/(4 * bat.pitch$N1))
for (k in 2:30){
  sure.vec[k] <- sure.grouplinear(bat.pitch$X1,1/(4 * bat.pitch$N1),nbreak=k)
}
khat.sure <- which.min(sure.vec) #2
delta.gl.sure <- if(khat.sure>1) grouplinear( bat.pitch$X1,1/(4 * bat.pitch$N1),nbreak=khat.sure) else spher(bat.pitch$X1,1/(4 * bat.pitch$N1))
tse.hat.delta.gl.sure <- sum(   (  ( bat.pitch$X2 - delta.gl.sure )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
tse.hat.delta.gl.sure/tse.hat.zero  # 0.174


## nonpitchers only

bat.npitch <- bat[bat$Pitcher.==0,]

# estimating TSE for various estimators
# run: functions.R(current folder), functions_XKB.R

ind.npitch <- bat.npitch$N2>10  # indicator for records with N2>=11 (among those with N1>=11)

tse.hat.zero <- sum(   (  ( bat.npitch$X2 - bat.npitch$X1 )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )

# grand mean
tse.hat.delta.gm <- sum(   (  ( bat.npitch$X2 - mean(bat.npitch$X1) )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
tse.hat.delta.gm/tse.hat.zero

# XKB theta.hat.M
delta.M <- thetahat.M(bat.npitch$X1,1/(4 * bat.npitch$N1))
tse.hat.delta.M <- sum(   (  ( bat.npitch$X2 - delta.M )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
tse.hat.delta.M/tse.hat.zero

# XKB theta.hat.SG
delta.SG <- thetahat.SG(bat.npitch$X1,1/(4 * bat.npitch$N1))
tse.hat.delta.SG <- sum(   (  ( bat.npitch$X2 - delta.SG )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
tse.hat.delta.SG/tse.hat.zero

# group-linear 
  # n^-1/3 bins
delta.gl <- grouplinear(x=bat.npitch$X1, v=1/(4 * bat.npitch$N1))
tse.hat.delta.gl <- sum(   (  ( bat.npitch$X2 - delta.gl )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
tse.hat.delta.gl/tse.hat.zero
  # oracle
rel.tse.breaks <- rep(NA,20)
delta.gl <- spher(x=bat.npitch$X1, v=1/(4 * bat.npitch$N1))
tse.hat.delta.gl <- sum(   (  ( bat.npitch$X2 - delta.gl )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
rel.tse.breaks[1] <- tse.hat.delta.gl/tse.hat.zero
for(i in 2:20){
  delta.gl <- grouplinear(x=bat.npitch$X1, v=1/(4 * bat.npitch$N1),nbreak = i)
  tse.hat.delta.gl <- sum(   (  ( bat.npitch$X2 - delta.gl )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
  rel.tse.breaks[i] <- tse.hat.delta.gl/tse.hat.zero
}

tse.hat.delta.gl.ol <- min(rel.tse.breaks) #.275
k.ol <- which.min(rel.tse.breaks) #1


# URE
# i) split into k intervals of equal length on log(v)
min.diff <- min(diff( sort(log( 1/(4 * bat.npitch$N1) )) )[diff( sort(log( 1/(4 * bat.npitch$N1) )) )>0])  # min_{i,j: v_i != v_j} |v_i-v_j|
kmax <- ceiling( diff(range(log( 1/(4 * bat.npitch$N1) )))/min.diff )
sure.vec <- rep(NA,kmax)

sure.vec[1] <- sure.spher(bat.npitch$X1,1/(4 * bat.npitch$N1))
for (k in 2:30){
  sure.vec[k] <- sure.grouplinear(bat.npitch$X1,1/(4 * bat.npitch$N1),nbreak=k)
}
khat.sure <- which.min(sure.vec) #2
delta.gl.sure <- if (khat.sure>1) grouplinear( bat.npitch$X1,1/(4 * bat.npitch$N1),nbreak=khat.sure) else spher(bat.npitch$X1,1/(4 * bat.npitch$N1))
tse.hat.delta.gl.sure <- sum(   (  ( bat.npitch$X2 - delta.gl.sure )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
tse.hat.delta.gl.sure/tse.hat.zero  # 0.338


## analog of XKB Figure 2

x <- bat$X1
v <- 1/(4 * bat$N1)

# XKB theta.hat.M
lambda.sure <- ifelse( g(0,X=x,A=v)*g(max(v)*1000,X=x,A=v) < 0, uniroot(g,c(0,max(v)*1000),X=x,A=v, tol=1e-9)$root, optim(c(mean(   pmax(  ( x-mean(x) )^2 - 	v,0  )   ), mean(x)),f,X=x,A=v, method = "L-BFGS-B",lower=c(0,-Inf))$par[1] )
mu.sure <- sum( v^2/(lambda.sure+v)^2 * x ) / sum( v^2/(lambda.sure+v)^2 )

# XKB theta.hat.SG
p <- length(x)
fit <- gpava( z = v, y = v * (1-1/p) / (x-mean(x))^2, weights = (x-mean(x))^2, solver = weighted.mean, ties="primary" )
bhat.sg <- pmin(  pmax( fit$x,0 ),1  )

# grouplinear
spher.bhat <- function(x.,v.){
n. <- length(x.)
if ( (n.==1) | (var(x.)==0) ) x. else {
	cstar <- max( 1-2*( max(v.)/mean(v.) )/(n.-1), 0)
	bhat <- min( cstar*mean(v.)/var(x.), 1 )
	rep(bhat,n.)  # returns a vector of equal elements
	}
}
spher.mean <- function(x.) rep(mean(x.),length(x.))

splitby <- cut(log(v),breaks=floor( n^(1/3) ), labels=F)
n <- length(x)
xsub <- split(x,splitby)
vsub <- split(v,splitby)
indexsub <- split(1:n,splitby)
bhatsub <- mapply(spher.bhat,xsub,vsub)
indexsub.unlist <- as.vector( unlist(indexsub) )
bhatsub.unlist <- as.vector( unlist(bhatsub) )
bhat <- bhatsub.unlist[order(indexsub.unlist)]
bhat.gl <- bhat

ahatsub <- lapply(xsub,spher.mean)
ahatsub.unlist <- as.vector( unlist(ahatsub) )
ahat <- ahatsub.unlist[order(indexsub.unlist)]

plot(bat$N1[order(bat$N1)], ( lambda.sure/(lambda.sure+v) )[order(bat$N1)], ylim=c(0,.9), xlab=expression(N[1][i] ~ " = Number of at-bats"),ylab='1 - (Shrinkage factor)',type='n',cex.lab=.65)
lines(bat$N1[order(bat$N1)], ( lambda.sure/(lambda.sure+v) )[order(bat$N1)], type='l', lty=3)
# lines(bat$N1[order(bat$N1)], ( 1-bhat.sg )[order(bat$N1)],type='l', lty=1)
#monotonize sg:
lines(sort(bat$N1), sort( 1-bhat.sg ), type='l', lty=1)
# lines(bat$N1[order(bat$N1)], ( 1-bhat.gl )[order(bat$N1)], type='l', lty=1)
legend("topleft", inset=.04, c(expression("SURE " ~ hat(theta)^M), expression("SP SURE " ~ hat(theta)^SG), expression("Group-linear " ~ hat(theta)^GL)), lty=c(3,2,1),cex=.5, bty='n')

# # eps for Biometrika
#
# setEPS()
# postscript("whatever-033115.eps")
# #plot(rnorm(100), main="Hey Some Data")
# dev.off()



# abline(h=mu.sure, col='blue', lty=3)
# lines( bat$N1[order(bat$N1)], ahat[order(bat$N1)], col='blue' )



###################
#   August 2015   #
###################

#dataframe 'bat' includes all players (with >10 at-bats in 1st half-season)
#'bat.pitch' includes only pitchers; dataframe 'bat.npitch' only nonpitchers


## plot estimate vs. X2
plot(bat.npitch$X2[ind],delta.gl[ind], cex=.25, xlab=expression(X[2] ~ '= transformed batting avg (2nd half)'),ylab=expression("predicted value"))
points(bat.npitch$X2[ind],thetahat.SG(bat.npitch$X1,1/(4 * bat.npitch$N1))[ind],col='red',cex=.25)
points(bat.npitch$X2,bat.npitch$X2,col='darkgreen',cex=.25)
# abline(a=0,b=1,lty=2,lwd=1)

## 'Pitcher.' and sampling variance are correlated:
cor(1/(4*bat$N1),bat$Pitcher.)

## examine how prior VARIANCE is correlated with number of at-bats
x <- bat$X1
v <- 1/(4 * bat$N1)
f <- function(x.,v.) max(0,var(x.)-mean(v.)) # estimate prior variance (for use within blocks)
splitby <- cut(log(v),breaks=floor( n^(1/3) ), labels=F)
n <- length(x)
xsub <- split(x,splitby)
vsub <- split(v,splitby)
prior.vars <- as.vector( mapply(f,xsub,vsub) )
t.test(prior.vars)  # p-val = 0.06

  ## within nonpitchers
x <- bat.npitch$X1
v <- 1/(4 * bat.npitch$N1)
f <- function(x.,v.) max(0,var(x.)-mean(v.)) # estimate prior variance (for use within blocks)
splitby <- cut(log(v),breaks=floor( n^(1/3) ), labels=F)
n <- length(x)
xsub <- split(x,splitby)
vsub <- split(v,splitby)
prior.vars <- as.vector( mapply(f,xsub,vsub) )
t.test(prior.vars)  # p-val = 0.06; about the same..


# shuffling

# a function that permutes hits (only for players with N1>10 )
bat.perm <- function(){
  bat <- bat.raw
  bat$N1 <- bat$AB.4. + bat$AB.5. + bat$AB.6.  # total number at-bats for 1st period
  bat$N2 <- bat$AB.7. + bat$AB.8. + bat$AB.9.10.  # total number at-bats for 2nd period
  bat$H1 <- bat$H.4. + bat$H.5. + bat$H.6.  # total number hits for 1st period
  bat$H2 <- bat$H.7. + bat$H.8. + bat$H.9.10.  # total number hits for 2nd period
  # bat$R1 <- bat$H1/bat$N1  # batting avg for 1st period
  # bat$R2 <- bat$H2/bat$N2  # batting avg for 2nd period
  # bat$X1 <- asin(  sqrt( (bat$H1+1/4)/(bat$N1+1/2) )  )  # transformed batting avg for 1st period
  # bat$X2 <- asin(  sqrt( (bat$H2+1/4)/(bat$N2+1/2) )  )  # transformed batting avg for 2nd period
  bat <- bat[bat$N1 > 10,]  # keep only records with N1>=11  
  
  bat$H1.perm <- NA
  for(i in 1:dim(bat)[1]){
    bat$H1.perm[i] <- rhyper(nn=1,m=bat$H1[i] + bat$H2[i],n=bat$N1[i] + bat$N2[i] -bat$H1[i] - bat$H2[i],k=bat$N1[i])
  }
  bat$H2.perm <- bat$H1 + bat$H2 - bat$H1.perm
  # head(cbind(bat$H1,bat$H1.perm,bat$H2,bat$H2.perm))
  bat$H1 <- bat$H1.perm 
  bat$H2 <- bat$H2.perm 
  
  bat$R1 <- bat$H1/bat$N1  # batting avg for 1st period
  bat$R2 <- bat$H2/bat$N2  # batting avg for 2nd period
  bat$X1 <- asin(  sqrt( (bat$H1+1/4)/(bat$N1+1/2) )  )  # transformed batting avg for 1st period
  bat$X2 <- asin(  sqrt( (bat$H2+1/4)/(bat$N2+1/2) )  )  # transformed batting avg for 2nd period
  bat <-  bat[,c('First.Name','Last.Name','Pitcher.','N1','N2','H1','H2','X1','X2')]
}

N <- 1000# num shuffling rounds

## all batters

tse.gm <- rep(NA,N)
tse.M <- rep(NA,N)
tse.SG <- rep(NA,N)
tse.gl <- rep(NA,N)
tse.gl.ol <- rep(NA,N)
tse.gl.sure <- rep(NA,N)

for(j in 1:N){
  bat <- bat.perm()
  # estimating TSE for various estimators
  # run: functions.R(current folder), functions_XKB.R
  ind <- bat$N2>10  # indicator for records with N2>=11 (among those with N1>=11)
  tse.hat.zero <- sum(   (  ( bat$X2 - bat$X1 )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  
  # grand mean
  tse.hat.delta.gm <- sum(   (  ( bat$X2 - mean(bat$X1) )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gm[j] <- tse.hat.delta.gm/tse.hat.zero
  
  # XKB theta.hat.M
  delta.M <- thetahat.M(bat$X1,1/(4 * bat$N1))
  tse.hat.delta.M <- sum(   (  ( bat$X2 - delta.M )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.M[j] <- tse.hat.delta.M/tse.hat.zero
  
  # XKB theta.hat.SG
  delta.SG <- thetahat.SG(bat$X1,1/(4 * bat$N1))
  tse.hat.delta.SG <- sum(   (  ( bat$X2 - delta.SG )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.SG[j] <- tse.hat.delta.SG/tse.hat.zero
  
  # group-linear  
    # num bins = n^1/3
  delta.gl <- grouplinear(x=bat$X1, v=1/(4 * bat$N1))
  tse.hat.delta.gl <- sum(   (  ( bat$X2 - delta.gl )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl[j] <- tse.hat.delta.gl/tse.hat.zero
    # oracle
  rel.tse.breaks <- rep(NA,20)
  delta.gl <- spher(x=bat$X1, v=1/(4 * bat$N1))
  tse.hat.delta.gl <- sum(   (  ( bat$X2 - delta.gl )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  rel.tse.breaks[1] <- tse.hat.delta.gl/tse.hat.zero
  for(i in 2:20){
    delta.gl <- grouplinear(x=bat$X1, v=1/(4 * bat$N1),nbreak = i)
    tse.hat.delta.gl <- sum(   (  ( bat$X2 - delta.gl )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
    rel.tse.breaks[i] <- tse.hat.delta.gl/tse.hat.zero
  }
  tse.gl.ol[j] <- min(rel.tse.breaks)
#   k.ol <- which.min(rel.tse.breaks)  
  
    # URE
  # i) split into k intervals of equal length on log(v)
  min.diff <- min(diff( sort(log( 1/(4 * bat$N1) )) )[diff( sort(log( 1/(4 * bat$N1) )) )>0])  # min_{i,j: v_i != v_j} |v_i-v_j|
  kmax <- ceiling( diff(range(log( 1/(4 * bat$N1) )))/min.diff )
  sure.vec <- rep(NA,kmax)
  
  sure.vec[1] <- sure.spher(bat$X1,1/(4 * bat$N1))
  for (k in 2:30){
    sure.vec[k] <- sure.grouplinear(bat$X1,1/(4 * bat$N1),nbreak=k)
  }
  khat.sure <- which.min(sure.vec)
  delta.gl.sure <- if(khat.sure>1) grouplinear( bat$X1,1/(4 * bat$N1),nbreak=khat.sure) else spher( bat$X1,1/(4 * bat$N1))
  tse.hat.delta.gl.sure <- sum(   (  ( bat$X2 - delta.gl.sure )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.sure[j] <- tse.hat.delta.gl.sure/tse.hat.zero
print(j)
}

tse.gm.all <- mean(tse.gm)
tse.M.all <- mean(tse.M)
tse.SG.all <- mean(tse.SG)
tse.gl.all <- mean(tse.gl)
tse.gl.ol.all <- mean(tse.gl.ol)
tse.gl.sure.all <- mean(tse.gl.sure)


## pitchers only

tse.gm <- rep(NA,N)
tse.M <- rep(NA,N)
tse.SG <- rep(NA,N)
tse.gl <- rep(NA,N)
tse.gl.ol <- rep(NA,N)
tse.gl.sure <- rep(NA,N)

for(j in 1:N){
  bat <- bat.perm()
  bat.pitch <- bat[bat$Pitcher.==1,]
  # estimating TSE for various estimators
  # run: functions.R(current folder), functions_XKB.R
  ind.pitch <- bat.pitch$N2>10  # indicator for records with N2>=11 (among those with N1>=11)
  tse.hat.zero <- sum(   (  ( bat.pitch$X2 - bat.pitch$X1 )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
  
  # grand mean
  tse.hat.delta.gm <- sum(   (  ( bat.pitch$X2 - mean(bat.pitch$X1) )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
  tse.gm[j] <- tse.hat.delta.gm/tse.hat.zero
  
  # XKB theta.hat.M
  delta.M <- thetahat.M(bat.pitch$X1,1/(4 * bat.pitch$N1))
  tse.hat.delta.M <- sum(   (  ( bat.pitch$X2 - delta.M )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
  tse.M[j] <- tse.hat.delta.M/tse.hat.zero
  
  # XKB theta.hat.SG
  delta.SG <- thetahat.SG(bat.pitch$X1,1/(4 * bat.pitch$N1))
  tse.hat.delta.SG <- sum(   (  ( bat.pitch$X2 - delta.SG )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
  tse.SG[j] <- tse.hat.delta.SG/tse.hat.zero
  
  # group-linear  
    # num bins = n^1/3
  delta.gl <- grouplinear(x=bat.pitch$X1, v=1/(4 * bat.pitch$N1))
  tse.hat.delta.gl <- sum(   (  ( bat.pitch$X2 - delta.gl )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
  tse.gl[j] <- tse.hat.delta.gl/tse.hat.zero
    # oracle
  rel.tse.breaks <- rep(NA,20)
  delta.gl <- spher(x=bat.pitch$X1, v=1/(4 * bat.pitch$N1))
  tse.hat.delta.gl <- sum(   (  ( bat.pitch$X2 - delta.gl )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
  rel.tse.breaks[1] <- tse.hat.delta.gl/tse.hat.zero
  for(i in 2:20){
    delta.gl <- grouplinear(x=bat.pitch$X1, v=1/(4 * bat.pitch$N1),nbreak = i)
    tse.hat.delta.gl <- sum(   (  ( bat.pitch$X2 - delta.gl )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
    rel.tse.breaks[i] <- tse.hat.delta.gl/tse.hat.zero
  }
  tse.gl.ol[j] <- min(rel.tse.breaks)
  #   k.ol <- which.min(rel.tse.breaks)  
  
    # URE
  # i) split into k intervals of equal length on log(v)
  min.diff <- min(diff( sort(log( 1/(4 * bat.pitch$N1) )) )[diff( sort(log( 1/(4 * bat.pitch$N1) )) )>0])  # min_{i,j: v_i != v_j} |v_i-v_j|
  kmax <- ceiling( diff(range(log( 1/(4 * bat.pitch$N1) )))/min.diff )
  sure.vec <- rep(NA,kmax)
  
  sure.vec[1] <- sure.spher(bat.pitch$X1,1/(4 * bat.pitch$N1))
  for (k in 2:kmax){
    sure.vec[k] <- sure.grouplinear(bat.pitch$X1,1/(4 * bat.pitch$N1),nbreak=k)
  }
  khat.sure <- which.min(sure.vec)
  delta.gl.sure <- if(khat.sure>1) grouplinear( bat.pitch$X1,1/(4 * bat.pitch$N1),nbreak=khat.sure) else spher(bat.pitch$X1,1/(4 * bat.pitch$N1))
  tse.hat.delta.gl.sure <- sum(   (  ( bat.pitch$X2 - delta.gl.sure )^2 - 1/ ( 4 * bat.pitch$N2 )  )[ind.pitch]   )
  tse.gl.sure[j] <- tse.hat.delta.gl.sure/tse.hat.zero
  print(j)  
}

(tse.gm.pitch <- mean(tse.gm))
(tse.M.pitch <- mean(tse.M))
(tse.SG.pitch <- mean(tse.SG))
(tse.gl.pitch <- mean(tse.gl))
(tse.gl.ol.pitch <- mean(tse.gl.ol))
(tse.gl.sure.pitch <- mean(tse.gl.sure))


## non-pitchers only

tse.gm <- rep(NA,N)
tse.M <- rep(NA,N)
tse.SG <- rep(NA,N)
tse.gl <- rep(NA,N)
tse.gl.sure <- rep(NA,N)

for(j in 1:N){
  bat <- bat.perm()
  bat.npitch <- bat[bat$Pitcher.==0,]
  # estimating TSE for various estimators
  # run: functions.R(current folder), functions_XKB.R
  ind.npitch <- bat.npitch$N2>10  # indicator for records with N2>=11 (among those with N1>=11)
  tse.hat.zero <- sum(   (  ( bat.npitch$X2 - bat.npitch$X1 )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
  
  # grand mean
  tse.hat.delta.gm <- sum(   (  ( bat.npitch$X2 - mean(bat.npitch$X1) )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
  tse.gm[j] <- tse.hat.delta.gm/tse.hat.zero
  
  # XKB theta.hat.M
  delta.M <- thetahat.M(bat.npitch$X1,1/(4 * bat.npitch$N1))
  tse.hat.delta.M <- sum(   (  ( bat.npitch$X2 - delta.M )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
  tse.M[j] <- tse.hat.delta.M/tse.hat.zero
  
  # XKB theta.hat.SG
  delta.SG <- thetahat.SG(bat.npitch$X1,1/(4 * bat.npitch$N1))
  tse.hat.delta.SG <- sum(   (  ( bat.npitch$X2 - delta.SG )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
  tse.SG[j] <- tse.hat.delta.SG/tse.hat.zero
  
  # group-linear  
    # num bins = n^1/3
  delta.gl <- grouplinear(x=bat.npitch$X1, v=1/(4 * bat.npitch$N1))
  tse.hat.delta.gl <- sum(   (  ( bat.npitch$X2 - delta.gl )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
  tse.gl[j] <- tse.hat.delta.gl/tse.hat.zero
  # oracle
  rel.tse.breaks <- rep(NA,20)
  delta.gl <- spher(x=bat.npitch$X1, v=1/(4 * bat.npitch$N1))
  tse.hat.delta.gl <- sum(   (  ( bat.npitch$X2 - delta.gl )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
  rel.tse.breaks[1] <- tse.hat.delta.gl/tse.hat.zero
  for(i in 2:20){
    delta.gl <- grouplinear(x=bat.npitch$X1, v=1/(4 * bat.npitch$N1),nbreak = i)
    tse.hat.delta.gl <- sum(   (  ( bat.npitch$X2 - delta.gl )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
    rel.tse.breaks[i] <- tse.hat.delta.gl/tse.hat.zero
  }
  tse.gl.ol[j] <- min(rel.tse.breaks)
  #   k.ol <- which.min(rel.tse.breaks)  
  
  # URE
  # i) split into k intervals of equal length on log(v)
  min.diff <- min(diff( sort(log( 1/(4 * bat.npitch$N1) )) )[diff( sort(log( 1/(4 * bat.npitch$N1) )) )>0])  # min_{i,j: v_i != v_j} |v_i-v_j|
  kmax <- ceiling( diff(range(log( 1/(4 * bat.npitch$N1) )))/min.diff )
  sure.vec <- rep(NA,kmax)
  
  sure.vec[1] <- sure.spher(bat.npitch$X1,1/(4 * bat.npitch$N1))
  for (k in 2:30){
    sure.vec[k] <- sure.grouplinear(bat.npitch$X1,1/(4 * bat.npitch$N1),nbreak=k)
  }
  khat.sure <- which.min(sure.vec)
  delta.gl.sure <- if(khat.sure>1) grouplinear( bat.npitch$X1,1/(4 * bat.npitch$N1),nbreak=khat.sure) else spher(bat.npitch$X1,1/(4 * bat.npitch$N1))
  tse.hat.delta.gl.sure <- sum(   (  ( bat.npitch$X2 - delta.gl.sure )^2 - 1/ ( 4 * bat.npitch$N2 )  )[ind.npitch]   )
  tse.gl.sure[j] <- tse.hat.delta.gl.sure/tse.hat.zero
print(j)
}

tse.gm.npitch <- mean(tse.gm)
tse.M.npitch <- mean(tse.M)
tse.SG.npitch <- mean(tse.SG)
tse.gl.npitch <- mean(tse.gl)
tse.gl.ol.npitch <- mean(tse.gl.ol)
tse.gl.sure.npitch <- mean(tse.gl.sure)

# summary
summary <- rbind(
gm=c(tse.gm.all,tse.gm.pitch,tse.gm.npitch),
M=c(tse.M.all,tse.M.pitch,tse.M.npitch),
SG=c(tse.SG.all,tse.SG.pitch,tse.SG.npitch),
gl=c(tse.gl.all,tse.gl.pitch,tse.gl.npitch),
gl.ol=c(tse.gl.ol.all,tse.gl.ol.pitch,tse.gl.ol.npitch),
gl.sure=c(tse.gl.sure.all,tse.gl.sure.pitch,tse.gl.sure.npitch)
)
( summary <- data.frame(all=summary[,1],pitch=summary[,2],npitch=summary[,3],row.names = c('gm','M','SG','gl','gl.ol','gl.sure')) )




