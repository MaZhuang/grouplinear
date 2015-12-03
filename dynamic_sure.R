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
			b[i,j]=sure
			for (k in i:(j-1)){
				#temp=sure.spher(x[i:k], v[i:k])+sure.spher(x[(k+1):j], v[(k+1):j])
				temp=b[i,k]+b[k+1,j]
				if (b[i,j]>temp){
				#	print(sure)
					a[i,j]=k
					b[i,j]=temp
				}
			}
			
		}
	}
	list(a,b)
}

DynamicSureMin=function(x,v,d=40){
    n <- length(x)
    a=matrix(rep(0,n*n),ncol=n) ##separation
	b=a ##value
	for (i in 1:(n-d+1)){
		j=i+d-1
		a[i,j]=j
		b[i,j]=sure.spher(x[i:j], v[i:j])
	}
	for (l in d:(n-1)){
		if (l %% 20==0){
		    print(l)
		}
		for (i in 1:(n-l)){
			j=l+i
			sure=sure.spher(x[i:j], v[i:j])
		#	print(sure)
			a[i,j]=j
			b[i,j]=sure
			if ((i+d-1)<=(j-d)){
			for (k in (i+d-1):(j-d)){
				#temp=sure.spher(x[i:k], v[i:k])+sure.spher(x[(k+1):j], v[(k+1):j])
				temp=b[i,k]+b[k+1,j]
				if (b[i,j]>temp){
				#	print(sure)
					a[i,j]=k
					b[i,j]=temp
				}
			}}
			
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

GroupSureMin<- function(x,v,d){ 
   c=DynamicSureMin(x,v,d)
   position=c[[1]]
   n=dim(position)[1]
   group=partition(position,1,n)
   group=c(0, group,n)
   group=unique(group)
   est=dynamic.grouplinear(x,v,group)
   return(est)
}
