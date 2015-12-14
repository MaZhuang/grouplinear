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

N <- 1# num shuffling rounds
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
  
  # James-Stein
  delta.JS <- JS(bat$X1,1/(4 * bat$N1))
  tse.hat.delta.JS <- sum(   (  ( bat$X2 - delta.JS )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
  tse.JS[j] <- tse.hat.delta.JS/tse.hat.zero
  
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
  
  
    
  # dynamic_group_linear_all_division
delta.dynamic=GroupSure(bat$X1,1/(4 * bat$N1))
tse.hat.delta.dynamic <- sum(   (  ( bat$X2 - delta.dynamic )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.gl.dynamic[j] =tse.hat.delta.dynamic/tse.hat.zero

  # dynamic_group_linear_with_minimum_bin_size_constraint
delta.dynamicMin=GroupSureMin(bat$X1,1/(4 * bat$N1),n^(2/3)*0.8)
tse.hat.delta.dynamicMin <- sum(   (  ( bat$X2 - delta.dynamicMin )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.gl.dynamicMin[j] =tse.hat.delta.dynamicMin/tse.hat.zero

delta.dynamicMin2=GroupSureMin(bat$X1,1/(4 * bat$N1),n^(2/3))
tse.hat.delta.dynamicMin2 <- sum(   (  ( bat$X2 - delta.dynamicMin2 )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.gl.dynamicMin2[j] =tse.hat.delta.dynamicMin2/tse.hat.zero

delta.dynamicMin3=GroupSureMin(bat$X1,1/(4 * bat$N1),n^(2/3)*1.2)
tse.hat.delta.dynamicMin3 <- sum(   (  ( bat$X2 - delta.dynamicMin3 )^2 - 1/ ( 4 * bat$N2 )  )[ind]   )
tse.gl.dynamicMin3[j] =tse.hat.delta.dynamicMin3/tse.hat.zero



####################pitchers only
bat_p <- bat[bat$Pitcher.==1,]
n_p=dim(bat_p)[1]
ind_p <- bat_p$N2>10  # indicator for records with N2>=11 (among those with N1>=11)

  tse.hat.zero_p <- sum(   (  ( bat_p$X2 - bat_p$X1 )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  
  # grand mean
  tse.hat.delta.gm_p <- sum(   (  ( bat_p$X2 - mean(bat_p$X1) )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.gm_p[j] <- tse.hat.delta.gm_p/tse.hat.zero_p
  
  # James-Stein
  delta.JS_p <- JS(bat_p$X1,1/(4 * bat_p$N1))
  tse.hat.delta.JS_p <- sum(   (  ( bat_p$X2 - delta.JS_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.JS_p[j] <- tse.hat.delta.JS_p/tse.hat.zero_p
  
  # XKB theta.hat.M
  delta.M_p <- thetahat.M(bat_p$X1,1/(4 * bat_p$N1))
  tse.hat.delta.M_p <- sum(   (  ( bat_p$X2 - delta.M_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.M_p[j] <- tse.hat.delta.M_p/tse.hat.zero_p
  
  # XKB theta.hat.SG
  delta.SG_p <- thetahat.SG(bat_p$X1,1/(4 * bat_p$N1))
  tse.hat.delta.SG_p <- sum(   (  ( bat_p$X2 - delta.SG_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.SG_p[j] <- tse.hat.delta.SG_p/tse.hat.zero_p
  
  # group-linear  
    # num bins = n^1/3
  delta.gl_p <- grouplinear(x=bat_p$X1, v=1/(4 * bat_p$N1))
  tse.hat.delta.gl_p <- sum(   (  ( bat_p$X2 - delta.gl_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.gl_p[j] <- tse.hat.delta.gl_p/tse.hat.zero_p
   
   
       # URE
  # i) split into k intervals of equal length on log(v)
  min.diff_p <- min(diff( sort(log( 1/(4 * bat_p$N1) )) )[diff( sort(log( 1/(4 * bat_p$N1) )) )>0])  # min_{i,j: v_i != v_j} |v_i-v_j|
  kmax_p <- ceiling( diff(range(log( 1/(4 * bat_p$N1) )))/min.diff_p )
  sure.vec_p <- rep(NA,kmax_p)
  
  sure.vec_p[1] <- sure.spher(bat_p$X1,1/(4 * bat_p$N1))
  for (k in 2:30){
    sure.vec_p[k] <- sure.grouplinear(bat_p$X1,1/(4 * bat_p$N1),nbreak=k)
  }
  khat.sure_p <- which.min(sure.vec_p)
  delta.gl.sure_p <- if(khat.sure_p>1) grouplinear( bat_p$X1,1/(4 * bat_p$N1),nbreak=khat.sure_p) else spher( bat_p$X1,1/(4 * bat_p$N1))
  tse.hat.delta.gl.sure_p <- sum(   (  ( bat_p$X2 - delta.gl.sure_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
  tse.gl.sure_p[j] <- tse.hat.delta.gl.sure_p/tse.hat.zero_p
  
  
  
   
  # dynamic_group_linear_all_division
delta.dynamic_p=GroupSure(bat_p$X1,1/(4 * bat_p$N1))
tse.hat.delta.dynamic_p <- sum(   (  ( bat_p$X2 - delta.dynamic_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
tse.gl.dynamic_p[j] =tse.hat.delta.dynamic_p/tse.hat.zero_p

  # dynamic_group_linear_with_minimum_bin_size_constraint
delta.dynamicMin_p=GroupSureMin(bat_p$X1,1/(4 * bat_p$N1),n_p^(2/3)*0.8)
tse.hat.delta.dynamicMin_p <- sum(   (  ( bat_p$X2 - delta.dynamicMin_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
tse.gl.dynamicMin_p[j] =tse.hat.delta.dynamicMin_p/tse.hat.zero_p

delta.dynamicMin2_p=GroupSureMin(bat_p$X1,1/(4 * bat_p$N1),n_p^(2/3))
tse.hat.delta.dynamicMin2_p <- sum(   (  ( bat_p$X2 - delta.dynamicMin2_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
tse.gl.dynamicMin2_p[j] =tse.hat.delta.dynamicMin2_p/tse.hat.zero_p

delta.dynamicMin3_p=GroupSureMin(bat_p$X1,1/(4 * bat_p$N1),n_p^(2/3)*1.2)
tse.hat.delta.dynamicMin3_p <- sum(   (  ( bat_p$X2 - delta.dynamicMin3_p )^2 - 1/ ( 4 * bat_p$N2 )  )[ind_p]   )
tse.gl.dynamicMin3_p[j] =tse.hat.delta.dynamicMin3_p/tse.hat.zero_p



#########################Nonpitchers only
bat_n <- bat[bat$Pitcher.==0,]
n_n=dim(bat_n)[1]
ind_n <- bat_n$N2>10  # indicator for records with N2>=11 (among those with N1>=11)

  tse.hat.zero_n <- sum(   (  ( bat_n$X2 - bat_n$X1 )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  
  # grand mean
  tse.hat.delta.gm_n <- sum(   (  ( bat_n$X2 - mean(bat_n$X1) )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.gm_n[j] <- tse.hat.delta.gm_n/tse.hat.zero_n
  
  # James-Stein
  delta.JS_n <- JS(bat_n$X1,1/(4 * bat_n$N1))
  tse.hat.delta.JS_n <- sum(   (  ( bat_n$X2 - delta.JS_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.JS_n[j] <- tse.hat.delta.JS_n/tse.hat.zero_n
  
  # XKB theta.hat.M
  delta.M_n <- thetahat.M(bat_n$X1,1/(4 * bat_n$N1))
  tse.hat.delta.M_n <- sum(   (  ( bat_n$X2 - delta.M_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.M_n[j] <- tse.hat.delta.M_n/tse.hat.zero_n
  
  # XKB theta.hat.SG
  delta.SG_n <- thetahat.SG(bat_n$X1,1/(4 * bat_n$N1))
  tse.hat.delta.SG_n <- sum(   (  ( bat_n$X2 - delta.SG_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.SG_n[j] <- tse.hat.delta.SG_n/tse.hat.zero_n
  
  # group-linear  
    # num bins = n^1/3
  delta.gl_n <- grouplinear(x=bat_n$X1, v=1/(4 * bat_n$N1))
  tse.hat.delta.gl_n <- sum(   (  ( bat_n$X2 - delta.gl_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.gl_n[j] <- tse.hat.delta.gl_n/tse.hat.zero_n
   
  
  
         # URE
  # i) split into k intervals of equal length on log(v)
  min.diff_n <- min(diff( sort(log( 1/(4 * bat_n$N1) )) )[diff( sort(log( 1/(4 * bat_n$N1) )) )>0])  # min_{i,j: v_i != v_j} |v_i-v_j|
  kmax_n <- ceiling( diff(range(log( 1/(4 * bat_n$N1) )))/min.diff_n )
  sure.vec_n <- rep(NA,kmax_n)
  
  sure.vec_n[1] <- sure.spher(bat_n$X1,1/(4 * bat_n$N1))
  for (k in 2:30){
    sure.vec_n[k] <- sure.grouplinear(bat_n$X1,1/(4 * bat_n$N1),nbreak=k)
  }
  khat.sure_n <- which.min(sure.vec_n)
  delta.gl.sure_n <- if(khat.sure_n>1) grouplinear( bat_n$X1,1/(4 * bat_n$N1),nbreak=khat.sure_n) else spher( bat_n$X1,1/(4 * bat_n$N1))
  tse.hat.delta.gl.sure_n <- sum(   (  ( bat_n$X2 - delta.gl.sure_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
  tse.gl.sure_n[j] <- tse.hat.delta.gl.sure_n/tse.hat.zero_n
  
  
  
  
  # dynamic_group_linear_all_division
delta.dynamic_n=GroupSure(bat_n$X1,1/(4 * bat_n$N1))
tse.hat.delta.dynamic_n <- sum(   (  ( bat_n$X2 - delta.dynamic_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
tse.gl.dynamic_n[j] =tse.hat.delta.dynamic_n/tse.hat.zero_n

  # dynamic_group_linear_with_minimum_bin_size_constraint
delta.dynamicMin_n=GroupSureMin(bat_n$X1,1/(4 * bat_n$N1),n_n^(2/3)*0.8)
tse.hat.delta.dynamicMin_n <- sum(   (  ( bat_n$X2 - delta.dynamicMin_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
tse.gl.dynamicMin_n[j] =tse.hat.delta.dynamicMin_n/tse.hat.zero_n

delta.dynamicMin2_n=GroupSureMin(bat_n$X1,1/(4 * bat_n$N1),n_n^(2/3))
tse.hat.delta.dynamicMin2_n <- sum(   (  ( bat_n$X2 - delta.dynamicMin2_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
tse.gl.dynamicMin2_n[j] =tse.hat.delta.dynamicMin2_n/tse.hat.zero_n

delta.dynamicMin3_n=GroupSureMin(bat_n$X1,1/(4 * bat_n$N1),n_n^(2/3)*1.2)
tse.hat.delta.dynamicMin3_n <- sum(   (  ( bat_n$X2 - delta.dynamicMin3_n )^2 - 1/ ( 4 * bat_n$N2 )  )[ind_n]   )
tse.gl.dynamicMin3_n[j] =tse.hat.delta.dynamicMin3_n/tse.hat.zero_n
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
tse.gm.all_p,tse.JS.all_p,tse.M.all_p,tse.SG.all_p,tse.gl.all_p,tse.gl.sure.all_p,tse.gl.dynamic.all_p,tse.gl.dynamicMin.all_p,tse.gl.dynamicMin2.all_p,tse.gl.dynamicMin3.all_p,
tse.gm.all_n,tse.JS.all_n,tse.M.all_n,tse.SG.all_n,tse.gl.all_n,tse.gl.sure.all_n,tse.gl.dynamic.all_n,tse.gl.dynamicMin.all_n,tse.gl.dynamicMin2.all_n,tse.gl.dynamicMin3.all_n
)

error=cbind(tse.gm,tse.JS,tse.M,tse.SG,tse.gl,tse.gl.sure,tse.gl.dynamic,tse.gl.dynamicMin,tse.gl.dynamicMin2,tse.gl.dynamicMin3,
tse.gm_p, tse.JS_p, tse.M_p, tse.SG_p, tse.gl_p, tse.gl.sure_p, tse.gl.dynamic_p, tse.gl.dynamicMin_p, tse.gl.dynamicMin2_p, tse.gl.dynamicMin3_p,
tse.gm_n,tse.JS_n,tse.M_n, tse.SG_n, tse.gl_n, tse.gl.sure_n, tse.gl.dynamic_n, tse.gl.dynamicMin_n, tse.gl.dynamicMin2_n, tse.gl.dynamicMin3_n
)


names=rownames(average)
write.table(average, "~/desktop/average.txt",sep="\t",row.names=names)
write.table(error, "~/desktop/error.txt",sep="\t",row.names=FALSE)

sd=sqrt(diag(cov(error)))
write.table(sd, "~/desktop/sd.txt",sep="\t",row.names=names)








