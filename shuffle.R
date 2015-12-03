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
tse.M <- rep(NA,N)
tse.SG <- rep(NA,N)
tse.gl <- rep(NA,N)
tse.gl.ol <- rep(NA,N)
tse.gl.sure <- rep(NA,N)
tse.gl.dynamic <- rep(NA,N)
tse.gl.dynamicMin <- rep(NA,N)
tse.gl.dynamicMin2 <- rep(NA,N)
tse.gl.dynamicMin3 <- rep(NA,N)

path=getwd()
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
  
  
delta.dynamic=GroupSure(bat$X1,1/(4 * bat$N1))
tse.hat.delta.dynamic <- sum(   (  ( bat$X2 - delta.dynamic )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.gl.dynamic[j] =tse.hat.delta.dynamic/tse.hat.zero


delta.dynamicMin=GroupSureMin(bat$X1,1/(4 * bat$N1),40)
tse.hat.delta.dynamicMin <- sum(   (  ( bat$X2 - delta.dynamicMin )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.gl.dynamicMin[j] =tse.hat.delta.dynamicMin/tse.hat.zero

delta.dynamicMin2=GroupSureMin(bat$X1,1/(4 * bat$N1),50)
tse.hat.delta.dynamicMin2 <- sum(   (  ( bat$X2 - delta.dynamicMin2 )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.gl.dynamicMin2[j] =tse.hat.delta.dynamicMin/tse.hat.zero

delta.dynamicMin3=GroupSureMin(bat$X1,1/(4 * bat$N1),60)
tse.hat.delta.dynamicMin3 <- sum(   (  ( bat$X2 - delta.dynamicMin3 )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.gl.dynamicMin3[j] =tse.hat.delta.dynamicMin/tse.hat.zero
}

tse.gm.all <- mean(tse.gm)
tse.M.all <- mean(tse.M)
tse.SG.all <- mean(tse.SG)
tse.gl.all <- mean(tse.gl)
tse.gl.ol.all <- mean(tse.gl.ol)
tse.gl.sure.all <- mean(tse.gl.sure)
tse.gl.dynamic.all <- mean(tse.gl.dynamic)
tse.gl.dynamicMin.all <- mean(tse.gl.dynamicMin)
tse.gl.dynamicMin2.all <- mean(tse.gl.dynamicMin2)
tse.gl.dynamicMin3.all <- mean(tse.gl.dynamicMin3)


average=c(tse.gm.all,tse.M.all,tse.SG.all,tse.gl.all,tse.gl.ol.all,tse.gl.sure.all,tse.gl.dynamic.all,tse.gl.dynamicMin.all,tse.gl.dynamicMin2.all,tse.gl.dynamicMin3.all )
error=cbind(tse.gm,tse.M,tse.SG,tse.gl,tse.gl.ol,tse.gl.sure,tse.gl.dynamic,tse.gl.dynamicMin,tse.gl.dynamicMin2,tse.gl.dynamicMin3)

names=colnames(error)
write.table(average, "~/desktop/average.txt",sep="\t",row.names=names)
write.table(error, "~/desktop/error.txt",sep="\t",row.names=FALSE)
