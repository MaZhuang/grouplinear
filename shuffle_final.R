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
if(!exists("foo", mode="function")) source("functions.R")
if(!exists("foo", mode="function")) source("functions_XKB.R")
if(!exists("foo", mode="function")) source("dynamic_sure.R")

####all batters
tse.gm <- rep(NA,N)
tse.JS <- rep(NA,N)
tse.M <- rep(NA,N)
tse.SG <- rep(NA,N)
tse.gl <- rep(NA,N)
tse.gl.sure <- rep(NA,N)
tse.gl.dynamic <- rep(NA,N)
tse.gl.dynamicMin <- rep(NA,N)
tse.gl.dynamicMin2 <- rep(NA,N)
tse.gl.dynamicMin3 <- rep(NA,N)

####covariate: additive
tse.gl.add <- rep(NA,N)
tse.gl.sure.add<- rep(NA,N)
tse.gl.dynamic.add <- rep(NA,N)
tse.gl.dynamicMin.add <- rep(NA,N)
tse.gl.dynamicMin2.add <- rep(NA,N)
tse.gl.dynamicMin3.add <- rep(NA,N)

####covariate: interaction
tse.gl.interaction <- rep(NA,N)
tse.gl.sure.interaction<- rep(NA,N)
tse.gl.dynamic.interaction <- rep(NA,N)
tse.gl.dynamicMin.interaction <- rep(NA,N)
tse.gl.dynamicMin2.interaction <- rep(NA,N)
tse.gl.dynamicMin3.interaction <- rep(NA,N)

####pitchers
tse.gm_p <- rep(NA,N)
tse.JS_p <- rep(NA,N)
tse.M_p <- rep(NA,N)
tse.SG_p <- rep(NA,N)
tse.gl_p <- rep(NA,N)
tse.gl.sure_p <- rep(NA,N)
tse.gl.dynamic_p <- rep(NA,N)
tse.gl.dynamicMin_p <- rep(NA,N)
tse.gl.dynamicMin2_p <- rep(NA,N)
tse.gl.dynamicMin3_p <- rep(NA,N)

####non_pitchers
tse.gm_n <- rep(NA,N)
tse.JS_n <- rep(NA,N)
tse.M_n <- rep(NA,N)
tse.SG_n <- rep(NA,N)
tse.gl_n <- rep(NA,N)
tse.gl.sure_n<- rep(NA,N)
tse.gl.dynamic_n <- rep(NA,N)
tse.gl.dynamicMin_n <- rep(NA,N)
tse.gl.dynamicMin2_n <- rep(NA,N)
tse.gl.dynamicMin3_n <- rep(NA,N)




path=getwd() #setwd("/Users/assafweinstein/Dropbox/Research/Cunhui/Code/grouplinear_Zhuang")
datapath=paste(path,'Brown_batting_data.txt',sep='/')
bat.raw <- read.table(datapath, header=TRUE, sep=",", quote="")

for(j in 1:N){
  bat <- bat.perm()
  cat('Loop:',j)
  index=order(bat$N1,decreasing=TRUE)
  
  bat=bat[index,]
  n=dim(bat)[1]
  # estimating TSE for various estimators
  # run: functions.R(current folder), functions_XKB.R
  ind <- bat$N2>10  # indicator for records with N2>=11 (among those with N1>=11)
  tse.zero <- sum(   (  ( bat$X2 - bat$X1 )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  
  # grand mean
  tse.delta.gm <- sum(   (  ( bat$X2 - mean(bat$X1) )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gm[j] <- tse.delta.gm/tse.zero
  
  # James-Stein
  delta.JS <- JS(bat$X1,1/(4 * bat$N1))
  tse.delta.JS <- sum(   (  ( bat$X2 - delta.JS )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.JS[j] <- tse.delta.JS/tse.zero
  
  # XKB theta.hat.M
  delta.M <- thetahat.M(bat$X1,1/(4 * bat$N1))
  tse.delta.M <- sum(   (  ( bat$X2 - delta.M )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.M[j] <- tse.delta.M/tse.zero
  
  # XKB theta.hat.SG
  delta.SG <- thetahat.SG(bat$X1,1/(4 * bat$N1))
  tse.delta.SG <- sum(   (  ( bat$X2 - delta.SG )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.SG[j] <- tse.delta.SG/tse.zero
  
  # group_linear: num bins = n^1/3
  delta.gl <- grouplinear(x=bat$X1, v=1/(4 * bat$N1))
  tse.delta.gl <- sum(   (  ( bat$X2 - delta.gl )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl[j] <- tse.delta.gl/tse.zero
  
  # group_linear: sure(equal-bins)
  delta.gl.sure <- grouplinear.sure(bat$X1, 1/(4 * bat$N1), kmax=min(ceiling(n^(1/3)/0.8),n))
  tse.delta.gl.sure <- sum(   (  ( bat$X2 - delta.gl.sure )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.sure[j] <- tse.delta.gl.sure/tse.zero
  
  # dynamic_group_linear_all_division
  delta.dynamic=GroupSure(bat$X1,1/(4 * bat$N1))
  tse.delta.dynamic <- sum(   (  ( bat$X2 - delta.dynamic )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamic[j] =tse.delta.dynamic/tse.zero
  
  # dynamic_group_linear_with_minimum_bin_size_constraint
  delta.dynamicMin=GroupSureMin(bat$X1,1/(4 * bat$N1),n^(2/3)*0.8)
  tse.delta.dynamicMin <- sum(   (  ( bat$X2 - delta.dynamicMin )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamicMin[j] =tse.delta.dynamicMin/tse.zero
  
  delta.dynamicMin2=GroupSureMin(bat$X1,1/(4 * bat$N1),n^(2/3))
  tse.delta.dynamicMin2 <- sum(   (  ( bat$X2 - delta.dynamicMin2 )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamicMin2[j] =tse.delta.dynamicMin2/tse.zero
  
  delta.dynamicMin3=GroupSureMin(bat$X1,1/(4 * bat$N1),n^(2/3)*1.2)
  tse.delta.dynamicMin3 <- sum(   (  ( bat$X2 - delta.dynamicMin3 )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamicMin3[j] =tse.delta.dynamicMin3/tse.zero
  
  
  
  ####################covariate: additive (1 + AB + pitcher)
  # group-linear methods only
  lm1 <- lm(bat$X1 ~ as.factor(bat$Pitcher.) + bat$N1, weights = 4 * bat$N1)
  y <- resid(lm1)
  sig.hat <- summary(lm1)$sigma
  #compute covariance of the residual vector
  cov.y <- diag(1/(4 * bat$N1)) - model.matrix(lm1) %*% (vcov(lm1)/sig.hat^2) %*% t(model.matrix(lm1)) # Note: division by sig.hat^2 makes a difference!
  s <- svd(cov.y)
  #sort(s$d)[1:5]
  dim <- nobs(lm1)-3
  u1 <- (s$u)[,1:dim]
  # transform (to a vector with full-rank diagonal cov matrix)
  y.white <- crossprod(u1,y)
  estimate.white <- grouplinear.zero(x=y.white, v=s$d[1:dim]) #tse=.22 if grouplinear() instead of grouplinear.zero used
  estimate.white.sure <- grouplinear.sure.zero(x=y.white,v=s$d[1:dim],kmax=min(ceiling(dim^(1/3)/0.8),dim))
  estimate.white.dynamic <- GroupSure.zero(x=y.white,v=s$d[1:dim])
  estimate.white.dynamicMin <- GroupSureMin.zero(x=y.white,v=s$d[1:dim],n^(2/3)*0.8)
  estimate.white.dynamicMin2 <- GroupSureMin.zero(x=y.white,v=s$d[1:dim],n^(2/3))
  estimate.white.dynamicMin3 <- GroupSureMin.zero(x=y.white,v=s$d[1:dim],n^(2/3)*1.2)
  # back-transform
  estimate <- u1 %*% estimate.white
  estimate.sure <- u1 %*% estimate.white.sure
  estimate.dynamic <- u1 %*% estimate.white.dynamic
  estimate.dynamicMin <- u1 %*% estimate.white.dynamicMin
  estimate.dynamicMin2 <- u1 %*% estimate.white.dynamicMin2
  estimate.dynamicMin3 <- u1 %*% estimate.white.dynamicMin3
  # final estimate
  delta.gl.add <- fitted(lm1) + estimate
  tse.delta.gl.add <- sum(   (  ( bat$X2 - delta.gl.add )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.add[j] <- tse.delta.gl.add/tse.zero #tse=.20
  
  delta.gl.sure.add <- fitted(lm1) + estimate.sure
  tse.delta.gl.sure.add <- sum(   (  ( bat$X2 - delta.gl.sure.add )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.sure.add[j] <- tse.delta.gl.sure.add/tse.zero #tse=.20
  
  delta.gl.dynamic.add <- fitted(lm1) + estimate.dynamic
  tse.delta.gl.dynamic.add <- sum(   (  ( bat$X2 - delta.gl.dynamic.add )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamic.add[j] <- tse.delta.gl.dynamic.add/tse.zero #tse=.20
  
  delta.gl.dynamicMin.add <- fitted(lm1) + estimate.dynamicMin
  tse.delta.gl.dynamicMin.add <- sum(   (  ( bat$X2 - delta.gl.dynamicMin.add )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamicMin.add[j] <- tse.delta.gl.dynamicMin.add/tse.zero #tse=.20
  
  delta.gl.dynamicMin2.add <- fitted(lm1) + estimate.dynamicMin2
  tse.delta.gl.dynamicMin2.add <- sum(   (  ( bat$X2 - delta.gl.dynamicMin2.add )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamicMin2.add[j] <- tse.delta.gl.dynamicMin2.add/tse.zero #tse=.20
  
  delta.gl.dynamicMin3.add <- fitted(lm1) + estimate.dynamicMin3
  tse.delta.gl.dynamicMin3.add <- sum(   (  ( bat$X2 - delta.gl.dynamicMin3.add )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamicMin3.add[j] <- tse.delta.gl.dynamicMin3.add/tse.zero #tse=.20
  
  
  
  ####################covariate: interaction (1 + AB*pitcher)
  # group-linear methods only
  lm1 <- lm(bat$X1 ~ as.factor(bat$Pitcher.) * bat$N1, weights = 4 * bat$N1)
  y <- resid(lm1)
  sig.hat <- summary(lm1)$sigma
  #compute covariance of the residual vector
  cov.y <- diag(1/(4 * bat$N1)) - model.matrix(lm1) %*% (vcov(lm1)/sig.hat^2) %*% t(model.matrix(lm1)) # Note: division by sig.hat^2 makes a difference!
  s <- svd(cov.y)
  #sort(s$d)[1:5]
  dim <- nobs(lm1)-4
  u1 <- (s$u)[,1:dim]
  # transform (to a vector with full-rank diagonal cov matrix)
  y.white <- crossprod(u1,y)
  estimate.white <- grouplinear.zero(x=y.white, v=s$d[1:dim]) #tse=.22 if grouplinear() instead of grouplinear.zero used
  estimate.white.sure <- grouplinear.sure.zero(x=y.white,v=s$d[1:dim],kmax=min(ceiling(dim^(1/3)/0.8),dim))
  estimate.white.dynamic <- GroupSure.zero(x=y.white,v=s$d[1:dim])
  estimate.white.dynamicMin <- GroupSureMin.zero(x=y.white,v=s$d[1:dim],n^(2/3)*0.8)
  estimate.white.dynamicMin2 <- GroupSureMin.zero(x=y.white,v=s$d[1:dim],d=n^(2/3)) #
  estimate.white.dynamicMin3 <- GroupSureMin.zero(x=y.white,v=s$d[1:dim],n^(2/3)*1.2)
  # back-transform
  estimate <- u1 %*% estimate.white
  estimate.sure <- u1 %*% estimate.white.sure
  estimate.dynamic <- u1 %*% estimate.white.dynamic
  estimate.dynamicMin <- u1 %*% estimate.white.dynamicMin
  estimate.dynamicMin2 <- u1 %*% estimate.white.dynamicMin2
  estimate.dynamicMin3 <- u1 %*% estimate.white.dynamicMin3
  
  # final estimate
  delta.gl.interaction <- fitted(lm1) + estimate
  tse.delta.gl.interaction <- sum(   (  ( bat$X2 - delta.gl.interaction )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.interaction[j] <- tse.delta.gl.interaction/tse.zero #tse=.20
  
  delta.gl.sure.interaction <- fitted(lm1) + estimate.sure
  tse.delta.gl.sure.interaction <- sum(   (  ( bat$X2 - delta.gl.sure.interaction )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.sure.interaction[j] <- tse.delta.gl.sure.interaction/tse.zero #tse=.20  
  
  delta.gl.dynamic.interaction <- fitted(lm1) + estimate.dynamic
  tse.delta.gl.dynamic.interaction <- sum(   (  ( bat$X2 - delta.gl.dynamic.interaction )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamic.interaction[j] <- tse.delta.gl.dynamic.interaction/tse.zero #tse=.20
  
  delta.gl.dynamicMin.interaction <- fitted(lm1) + estimate.dynamicMin
  tse.delta.gl.dynamicMin.interaction <- sum(   (  ( bat$X2 - delta.gl.dynamicMin.interaction )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamicMin.interaction[j] <- tse.delta.gl.dynamicMin.interaction/tse.zero #tse=.20
  
  delta.gl.dynamicMin2.interaction <- fitted(lm1) + estimate.dynamicMin2
  tse.delta.gl.dynamicMin2.interaction <- sum(   (  ( bat$X2 - delta.gl.dynamicMin2.interaction )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamicMin2.interaction[j] <- tse.delta.gl.dynamicMin2.interaction/tse.zero #tse=.20
  
  delta.gl.dynamicMin3.interaction <- fitted(lm1) + estimate.dynamicMin3
  tse.delta.gl.dynamicMin3.interaction <- sum(   (  ( bat$X2 - delta.gl.dynamicMin3.interaction )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamicMin3.interaction[j] <- tse.delta.gl.dynamicMin3.interaction/tse.zero #tse=.20
  
  
  
  
  ####################pitchers only
  bat_p <- bat[bat$Pitcher.==1,]
  n_p=dim(bat_p)[1]
  ind_p <- bat_p$N2>10  # indicator for records with N2>=11 (among those with N1>=11)
  
  tse.zero_p <- sum(   (  ( bat_p$X2 - bat_p$X1 )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  
  # grand mean
  tse.delta.gm_p <- sum(   (  ( bat_p$X2 - mean(bat_p$X1) )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.gm_p[j] <- tse.delta.gm_p/tse.zero_p
  
  # James-Stein
  delta.JS_p <- JS(bat_p$X1,1/(4 * bat_p$N1))
  tse.delta.JS_p <- sum(   (  ( bat_p$X2 - delta.JS_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.JS_p[j] <- tse.delta.JS_p/tse.zero_p
  
  # XKB theta.hat.M
  delta.M_p <- thetahat.M(bat_p$X1,1/(4 * bat_p$N1))
  tse.delta.M_p <- sum(   (  ( bat_p$X2 - delta.M_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.M_p[j] <- tse.delta.M_p/tse.zero_p
  
  # XKB theta.hat.SG
  delta.SG_p <- thetahat.SG(bat_p$X1,1/(4 * bat_p$N1))
  tse.delta.SG_p <- sum(   (  ( bat_p$X2 - delta.SG_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.SG_p[j] <- tse.delta.SG_p/tse.zero_p
  
  # group_linear: num bins = n^1/3
  delta.gl_p <- grouplinear(x=bat_p$X1, v=1/(4 * bat_p$N1))
  tse.delta.gl_p <- sum(   (  ( bat_p$X2 - delta.gl_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.gl_p[j] <- tse.delta.gl_p/tse.zero_p
  
  
  # group_linear: sure(equal-bins)
  delta.gl.sure_p <- grouplinear.sure(bat_p$X1, 1/(4 * bat_p$N1), kmax=min(ceiling(n_p^(1/3)/0.8),n_p))
  tse.delta.gl.sure_p <- sum(   (  ( bat_p$X2 - delta.gl.sure_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.gl.sure_p[j] <- tse.delta.gl.sure_p/tse.zero_p
    
  # dynamic_group_linear_all_division
  delta.dynamic_p=GroupSure(bat_p$X1,1/(4 * bat_p$N1))
  tse.delta.dynamic_p <- sum(   (  ( bat_p$X2 - delta.dynamic_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.gl.dynamic_p[j] =tse.delta.dynamic_p/tse.zero_p
  
  # dynamic_group_linear_with_minimum_bin_size_constraint
  delta.dynamicMin_p=GroupSureMin(bat_p$X1,1/(4 * bat_p$N1),n_p^(2/3)*0.8)
  tse.delta.dynamicMin_p <- sum(   (  ( bat_p$X2 - delta.dynamicMin_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.gl.dynamicMin_p[j] =tse.delta.dynamicMin_p/tse.zero_p
  
  delta.dynamicMin2_p=GroupSureMin(bat_p$X1,1/(4 * bat_p$N1),n_p^(2/3))
  tse.delta.dynamicMin2_p <- sum(   (  ( bat_p$X2 - delta.dynamicMin2_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.gl.dynamicMin2_p[j] =tse.delta.dynamicMin2_p/tse.zero_p
  
  delta.dynamicMin3_p=GroupSureMin(bat_p$X1,1/(4 * bat_p$N1),n_p^(2/3)*1.2)
  tse.delta.dynamicMin3_p <- sum(   (  ( bat_p$X2 - delta.dynamicMin3_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.gl.dynamicMin3_p[j] =tse.delta.dynamicMin3_p/tse.zero_p
  
  
  
  #########################Nonpitchers only
  bat_n <- bat[bat$Pitcher.==0,]
  n_n=dim(bat_n)[1]
  ind_n <- bat_n$N2>10  # indicator for records with N2>=11 (among those with N1>=11)
  
  tse.zero_n <- sum(   (  ( bat_n$X2 - bat_n$X1 )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  
  # grand mean
  tse.delta.gm_n <- sum(   (  ( bat_n$X2 - mean(bat_n$X1) )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.gm_n[j] <- tse.delta.gm_n/tse.zero_n
  
  # James-Stein
  delta.JS_n <- JS(bat_n$X1,1/(4 * bat_n$N1))
  tse.delta.JS_n <- sum(   (  ( bat_n$X2 - delta.JS_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.JS_n[j] <- tse.delta.JS_n/tse.zero_n
  
  # XKB theta.hat.M
  delta.M_n <- thetahat.M(bat_n$X1,1/(4 * bat_n$N1))
  tse.delta.M_n <- sum(   (  ( bat_n$X2 - delta.M_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.M_n[j] <- tse.delta.M_n/tse.zero_n
  
  # XKB theta.hat.SG
  delta.SG_n <- thetahat.SG(bat_n$X1,1/(4 * bat_n$N1))
  tse.delta.SG_n <- sum(   (  ( bat_n$X2 - delta.SG_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.SG_n[j] <- tse.delta.SG_n/tse.zero_n
  
  # group_linear: num bins = n^1/3
  delta.gl_n <- grouplinear(x=bat_n$X1, v=1/(4 * bat_n$N1))
  tse.delta.gl_n <- sum(   (  ( bat_n$X2 - delta.gl_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.gl_n[j] <- tse.delta.gl_n/tse.zero_n
    
  # group_linear: sure(equal-bins)
  delta.gl.sure_n <- grouplinear.sure(bat_n$X1,1/(4 * bat_n$N1),kmax=min(ceiling(n_n^(1/3)/0.8),n_n))
  tse.delta.gl.sure_n <- sum(   (  ( bat_n$X2 - delta.gl.sure_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.gl.sure_n[j] <- tse.delta.gl.sure_n/tse.zero_n
  
  # dynamic_group_linear_all_division
  delta.dynamic_n=GroupSure(bat_n$X1,1/(4 * bat_n$N1))
  tse.delta.dynamic_n <- sum(   (  ( bat_n$X2 - delta.dynamic_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.gl.dynamic_n[j] =tse.delta.dynamic_n/tse.zero_n
  
  # dynamic_group_linear_with_minimum_bin_size_constraint
  delta.dynamicMin_n=GroupSureMin(bat_n$X1,1/(4 * bat_n$N1),n_n^(2/3)*0.8)
  tse.delta.dynamicMin_n <- sum(   (  ( bat_n$X2 - delta.dynamicMin_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.gl.dynamicMin_n[j] =tse.delta.dynamicMin_n/tse.zero_n
  
  delta.dynamicMin2_n=GroupSureMin(bat_n$X1,1/(4 * bat_n$N1),n_n^(2/3))
  tse.delta.dynamicMin2_n <- sum(   (  ( bat_n$X2 - delta.dynamicMin2_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.gl.dynamicMin2_n[j] =tse.delta.dynamicMin2_n/tse.zero_n
  
  delta.dynamicMin3_n=GroupSureMin(bat_n$X1,1/(4 * bat_n$N1),n_n^(2/3)*1.2)
  tse.delta.dynamicMin3_n <- sum(   (  ( bat_n$X2 - delta.dynamicMin3_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.gl.dynamicMin3_n[j] =tse.delta.dynamicMin3_n/tse.zero_n
}

tse.gm.all <- mean(tse.gm)
tse.JS.all <- mean(tse.JS)
tse.M.all <- mean(tse.M)
tse.SG.all <- mean(tse.SG)
tse.gl.all <- mean(tse.gl)
tse.gl.sure.all <- mean(tse.gl.sure)
tse.gl.dynamic.all <- mean(tse.gl.dynamic)
tse.gl.dynamicMin.all <- mean(tse.gl.dynamicMin)
tse.gl.dynamicMin2.all <- mean(tse.gl.dynamicMin2)
tse.gl.dynamicMin3.all <- mean(tse.gl.dynamicMin3)

tse.gl.all.add <- mean(tse.gl.add)
tse.gl.sure.all.add <- tse.gl.sure.add
tse.gl.dynamic.all.add <- mean(tse.gl.dynamic.add)
tse.gl.dynamicMin.all.add <- mean(tse.gl.dynamicMin.add)
tse.gl.dynamicMin2.all.add <- mean(tse.gl.dynamicMin2.add)
tse.gl.dynamicMin3.all.add <- mean(tse.gl.dynamicMin3.add)

tse.gl.all.interaction <- mean(tse.gl.interaction)
tse.gl.sure.all.interaction <- tse.gl.sure.interaction
tse.gl.dynamic.all.interaction <- mean(tse.gl.dynamic.interaction)
tse.gl.dynamicMin.all.interaction <- mean(tse.gl.dynamicMin.interaction)
tse.gl.dynamicMin2.all.interaction <- mean(tse.gl.dynamicMin2.interaction)
tse.gl.dynamicMin3.all.interaction <- mean(tse.gl.dynamicMin3.interaction)

tse.gm.all_p <- mean(tse.gm_p)
tse.JS.all_p <- mean(tse.JS_p)
tse.M.all_p <- mean(tse.M_p)
tse.SG.all_p <- mean(tse.SG_p)
tse.gl.all_p <- mean(tse.gl_p)
tse.gl.sure.all_p <- mean(tse.gl.sure_p)
tse.gl.dynamic.all_p <- mean(tse.gl.dynamic_p)
tse.gl.dynamicMin.all_p <- mean(tse.gl.dynamicMin_p)
tse.gl.dynamicMin2.all_p <- mean(tse.gl.dynamicMin2_p)
tse.gl.dynamicMin3.all_p <- mean(tse.gl.dynamicMin3_p)


tse.gm.all_n <- mean(tse.gm_n)
tse.JS.all_n <- mean(tse.JS_n)
tse.M.all_n <- mean(tse.M_n)
tse.SG.all_n <- mean(tse.SG_n)
tse.gl.all_n <- mean(tse.gl_n)
tse.gl.sure.all_n <- mean(tse.gl.sure_n)
tse.gl.dynamic.all_n <- mean(tse.gl.dynamic_n)
tse.gl.dynamicMin.all_n <- mean(tse.gl.dynamicMin_n)
tse.gl.dynamicMin2.all_n <- mean(tse.gl.dynamicMin2_n)
tse.gl.dynamicMin3.all_n <- mean(tse.gl.dynamicMin3_n)



average=rbind(tse.gm.all,tse.JS.all,tse.M.all,tse.SG.all,tse.gl.all,tse.gl.sure.all,tse.gl.dynamic.all,tse.gl.dynamicMin.all,tse.gl.dynamicMin2.all,tse.gl.dynamicMin3.all,
              tse.gl.all.add,tse.gl.sure.all.add,tse.gl.dynamic.all.add,tse.gl.dynamicMin.all.add,tse.gl.dynamicMin2.all.add,tse.gl.dynamicMin3.all.add,
              tse.gl.all.interaction,tse.gl.sure.all.interaction,tse.gl.dynamic.all.interaction,tse.gl.dynamicMin.all.interaction,tse.gl.dynamicMin2.all.interaction,tse.gl.dynamicMin3.all.interaction,
              tse.gm.all_p,tse.JS.all_p,tse.M.all_p,tse.SG.all_p,tse.gl.all_p,tse.gl.sure.all_p,tse.gl.dynamic.all_p,tse.gl.dynamicMin.all_p,tse.gl.dynamicMin2.all_p,tse.gl.dynamicMin3.all_p,
              tse.gm.all_n,tse.JS.all_n,tse.M.all_n,tse.SG.all_n,tse.gl.all_n,tse.gl.sure.all_n,tse.gl.dynamic.all_n,tse.gl.dynamicMin.all_n,tse.gl.dynamicMin2.all_n,tse.gl.dynamicMin3.all_n
)

error=cbind(tse.gm,tse.JS,tse.M,tse.SG,tse.gl,tse.gl.sure,tse.gl.dynamic,tse.gl.dynamicMin,tse.gl.dynamicMin2,tse.gl.dynamicMin3,
            tse.gl.add,tse.gl.sure.add,tse.gl.dynamic.add,tse.gl.dynamicMin.add,tse.gl.dynamicMin2.add,tse.gl.dynamicMin3.add,
            tse.gl.interaction,tse.gl.sure.interaction,tse.gl.dynamic.interaction,tse.gl.dynamicMin.interaction,tse.gl.dynamicMin2.interaction,tse.gl.dynamicMin3.interaction,
            tse.gm_p, tse.JS_p, tse.M_p, tse.SG_p, tse.gl_p, tse.gl.sure_p, tse.gl.dynamic_p, tse.gl.dynamicMin_p, tse.gl.dynamicMin2_p, tse.gl.dynamicMin3_p,
            tse.gm_n,tse.JS_n,tse.M_n, tse.SG_n, tse.gl_n, tse.gl.sure_n, tse.gl.dynamic_n, tse.gl.dynamicMin_n, tse.gl.dynamicMin2_n, tse.gl.dynamicMin3_n
)
sd=sqrt(diag(cov(error)))



names=rownames(average)
average.path=paste(path,'average.txt',sep='/')
error.path=paste(path,'error.txt',sep='/')
sd.path=paste(path,'sd.txt',sep='/')



write.table(average, average.path, sep="\t",row.names=names)
write.table(error, error.path, sep="\t",row.names=FALSE)
write.table(sd, sd.path, sep="\t",row.names=names)



#write.table(average, "~/desktop/average.txt",sep="\t",row.names=names)
#write.table(error, "~/desktop/error.txt",sep="\t",row.names=FALSE)
#write.table(sd, "~/desktop/sd.txt",sep="\t",row.names=names)








