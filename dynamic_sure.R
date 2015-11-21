DynamicSure=function(x,v){
    n <- length(x)
    a=matrix(rep(0,n*n),ncol=n) ##separation
	b=a ##value
	for(i in 1:n){
		a[i,i]=i
		b[i,i]=v[i]
	}
	for (l in 1:(n-1)){
		if (l %% 20==0){
		    print(l)
		}
		for (i in 1:(n-l)){
			j=l+i
			sure=sure.spher(x[i:j], v[i:j])
		#	print(sure)
			a[i,j]=j
			for (k in i:(j-1)){
				temp=sure.spher(x[i:k], v[i:k])+sure.spher(x[(k+1):j], v[(k+1):j])
				if (sure>temp){
					sure=temp
				#	print(sure)
					a[i,j]=k
					b[i,j]=sure
				}
			}
			
		}
	}
	list(a,b)
}
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

partition=function(position,i,j){
	if (position[i,j]==j){
		return(j)
	}else if (position[i,j]==i){
		return(i)
	}
	else{
		#par=c(partition(position,i,position[i,j]),position[i,j])
		#par=c(par,partition(position,position[i,j],j))
		#return(par)
		a=partition(position,i,position[i,j])
		b=partition(position,position[i,j],j)
		return(c(a,position[i,j],b))
	}
}



dynamic.grouplinear <- function(x,v,group){ #nbreak=num of bins
	ngroup <- length(group)
	n=length(x)
	est=rep(0,n)
	for (i in 1:(ngroup-1)){
		est[(group[i]+1):group[i+1]]=spher(x[(group[i]+1):group[i+1]],v[(group[i]+1):group[i+1]])
	}
	est
}

