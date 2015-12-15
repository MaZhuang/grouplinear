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

N <- 3# num shuffling rounds
if(!exists("foo", mode="function")) source("functions.R")
if(!exists("foo", mode="function")) source("functions_XKB.R")
if(!exists("foo", mode="function")) source("dynamic_sure.R")

## all batters
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

tse.gl.add <- rep(NA,N)
tse.gl.dynamic.add <- rep(NA,N)
tse.gl.dynamicMin.add <- rep(NA,N)
tse.gl.dynamicMin2.add <- rep(NA,N)
tse.gl.dynamicMin3.add <- rep(NA,N)

tse.gl.interaction <- rep(NA,N)
tse.gl.dynamic.interaction <- rep(NA,N)
tse.gl.dynamicMin.interaction <- rep(NA,N)
tse.gl.dynamicMin2.interaction <- rep(NA,N)
tse.gl.dynamicMin3.interaction <- rep(NA,N)



path=getwd() # setwd("/Users/assafweinstein/Dropbox/Research/Cunhui/Code/grouplinear_Zhuang")
datapath=paste(path,'Brown_batting_data.txt',sep='/')
bat.raw <- read.table(datapath, header=TRUE, sep=",", quote="")

for(j in 1:N){
  bat <- bat.perm()
  cat('Loop:',j,' ')
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
  
  # group-linear  
  # num bins = n^1/3
  delta.gl <- grouplinear(x=bat$X1, v=1/(4 * bat$N1))
  tse.delta.gl <- sum(   (  ( bat$X2 - delta.gl )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl[j] <- tse.delta.gl/tse.zero
  
 
  # SURE(equal-bins)
  delta.gl.sure <- grouplinear.sure(x=bat$X1, v=1/(4 * bat$N1),kmax=60)
  tse.delta.gl.sure <- sum(   (  ( bat$X2 - delta.gl.sure )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.sure[j] <- tse.delta.gl.sure/tse.zero
  
  
  delta.dynamic=GroupSure(bat$X1,1/(4 * bat$N1))
  tse.delta.dynamic <- sum(   (  ( bat$X2 - delta.dynamic )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamic[j] =tse.delta.dynamic/tse.zero
  
  #position=c[[1]]
  #n=dim(position)[1]
  #group=partition(position,1,n)
  #group=c(0, group,n)
  #group=unique(group)  
  
  delta.dynamicMin=GroupSureMin(bat$X1,1/(4 * bat$N1),40)
  tse.delta.dynamicMin <- sum(   (  ( bat$X2 - delta.dynamicMin )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamicMin[j] =tse.delta.dynamicMin/tse.zero

  #cmin=DynamicSureMin(bat$X1,1/(4 * bat$N1),40)
  #position=cmin[[1]]
  #n=dim(position)[1]
  #group=partition(position,1,n)
  #group=c(0, group,n)
  #group=unique(group)
  
  delta.dynamicMin2=GroupSureMin(bat$X1,1/(4 * bat$N1),50)
  tse.delta.dynamicMin2 <- sum(   (  ( bat$X2 - delta.dynamicMin2 )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamicMin2[j] =tse.delta.dynamicMin2/tse.zero
  
  delta.dynamicMin3=GroupSureMin(bat$X1,1/(4 * bat$N1),60)
  tse.delta.dynamicMin3 <- sum(   (  ( bat$X2 - delta.dynamicMin3 )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.dynamicMin3[j] =tse.delta.dynamicMin3/tse.zero

  # Including covariates

  # 1 + AB + pitcher
    # group-linear
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
  estimate.white <- grouplinear.zero(x=y.white, v=s$d[1:564]) #tse=.22 if grouplinear() instead of grouplinear.zero used
  estimate.white.dynamic <- GroupSure.zero(x=y.white,v=s$d[1:564])
  estimate.white.dynamicMin <- GroupSureMin.zero(x=y.white,v=s$d[1:564],40)
  estimate.white.dynamicMin2 <- GroupSureMin.zero(x=y.white,v=s$d[1:564],50)
  estimate.white.dynamicMin3 <- GroupSureMin.zero(x=y.white,v=s$d[1:564],60)
  # back-transform
  estimate <- u1 %*% estimate.white
  estimate.dynamic <- u1 %*% estimate.white.dynamic
  estimate.dynamicMin <- u1 %*% estimate.white.dynamicMin
  estimate.dynamicMin2 <- u1 %*% estimate.white.dynamicMin2
  estimate.dynamicMin3 <- u1 %*% estimate.white.dynamicMin3
  # final estimate
  delta.gl.add <- fitted(lm1) + estimate
  tse.delta.gl.add <- sum(   (  ( bat$X2 - delta.gl.add )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.add[j] <- tse.delta.gl.add/tse.zero #tse=.20
  
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
  estimate.white <- grouplinear.zero(x=y.white, v=s$d[1:563]) #tse=.22 if grouplinear() instead of grouplinear.zero used
  estimate.white.dynamic <- GroupSure.zero(x=y.white,v=s$d[1:563])
  estimate.white.dynamicMin <- GroupSureMin.zero(x=y.white,v=s$d[1:563],40)
  estimate.white.dynamicMin2 <- GroupSureMin.zero(x=y.white,v=s$d[1:563],d=50) #
  estimate.white.dynamicMin3 <- GroupSureMin.zero(x=y.white,v=s$d[1:563],60)
    # back-transform
  estimate <- u1 %*% estimate.white
  estimate.dynamic <- u1 %*% estimate.white.dynamic
  estimate.dynamicMin <- u1 %*% estimate.white.dynamicMin
  estimate.dynamicMin2 <- u1 %*% estimate.white.dynamicMin2
  estimate.dynamicMin3 <- u1 %*% estimate.white.dynamicMin3
    # final estimate
  delta.gl.interaction <- fitted(lm1) + estimate
  tse.delta.gl.interaction <- sum(   (  ( bat$X2 - delta.gl.interaction )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.gl.interaction[j] <- tse.delta.gl.interaction/tse.zero #tse=.20
  
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
}

tse.gm.all <- mean(tse.gm)
tse.JS.all <- mean(tse.JS)
tse.M.all <- mean(tse.M)
tse.SG.all <- mean(tse.SG)
tse.gl.all <- mean(tse.gl)
tse.gl.dynamic.all <- mean(tse.gl.dynamic)
tse.gl.dynamicMin.all <- mean(tse.gl.dynamicMin)
tse.gl.dynamicMin2.all <- mean(tse.gl.dynamicMin2)
tse.gl.dynamicMin3.all <- mean(tse.gl.dynamicMin3)

tse.gl.add.all <- mean(tse.gl.add)
tse.gl.dynamic.add.all <- mean(tse.gl.dynamic.add)
tse.gl.dynamicMin.add.all <- mean(tse.gl.dynamicMin.add)
tse.gl.dynamicMin2.add.all <- mean(tse.gl.dynamicMin2.add)
tse.gl.dynamicMin3.add.all <- mean(tse.gl.dynamicMin3.add)

tse.gl.interaction.all <- mean(tse.gl.interaction)
tse.gl.dynamic.interaction.all <- mean(tse.gl.dynamic.interaction)
tse.gl.dynamicMin.interaction.all <- mean(tse.gl.dynamicMin.interaction)
tse.gl.dynamicMin2.interaction.all <- mean(tse.gl.dynamicMin2.interaction)
tse.gl.dynamicMin3.interaction.all <- mean(tse.gl.dynamicMin3.interaction)

average=c(tse.gm.all,tse.JS.all,tse.M.all,tse.SG.all,tse.gl.all,tse.gl.dynamic.all,tse.gl.dynamicMin.all,tse.gl.dynamicMin2.all,tse.gl.dynamicMin3.all,
          tse.gl.add.all,tse.gl.dynamic.add.all,tse.gl.dynamicMin.add.all,tse.gl.dynamicMin2.add.all,tse.gl.dynamicMin3.add.all,
          tse.gl.interaction.all,tse.gl.dynamic.interaction.all,tse.gl.dynamicMin.interaction.all,tse.gl.dynamicMin2.interaction.all,tse.gl.dynamicMin3.interaction.all
          )

error=cbind(tse.gm,tse.JS,tse.M,tse.SG,tse.gl,tse.gl.dynamic,tse.gl.dynamicMin,tse.gl.dynamicMin2,tse.gl.dynamicMin3,
            tse.gl.add,tse.gl.dynamic.add,tse.gl.dynamicMin.add,tse.gl.dynamicMin2.add,tse.gl.dynamicMin3.add,
            tse.gl.interaction,tse.gl.dynamic.interaction,tse.gl.dynamicMin.interaction,tse.gl.dynamicMin2.interaction,tse.gl.dynamicMin3.interaction
            )


names=colnames(error)
write.table(average, "~/desktop/average.txt",sep="\t",row.names=names)
write.table(error, "~/desktop/error.txt",sep="\t",row.names=FALSE)









