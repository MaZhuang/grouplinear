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



partition=function(position,i,j){
	if (position[i,j]==j){
		return(j)
	}else if (position[i,j]==i){
		return(i)
	}
	else{
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

GroupSure<- function(x,v){ 
   c=DynamicSure(x,v)
   position=c[[1]]
   n=dim(position)[1]
   group=partition(position,1,n)
   group=c(0, group,n)
   group=unique(group)
   est=dynamic.grouplinear(x,v,group)
   return(est)
}

