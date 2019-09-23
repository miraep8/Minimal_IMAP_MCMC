# This function takes in the lower triangular adjacency matrix
# and the permutation and returns the full adjacency matrix

adjpermtoDAG<-function(n,permu,lowertriadjy){
	DAGadjy<-mat.or.vec(n,n)	
	for (kk in 1:(n-1)){
		for (ll in (1+kk):n){
			if(lowertriadjy[ll,kk]==1){
				DAGadjy[permu[ll],permu[kk]]<-1
			}
		}
	}
	return(DAGadjy)
}

set.seed(333)

nodelabels<-sample(c(1:n))

lowertriadjy<-matrix(0,n,n)

lowertriadjy[lower.tri(lowertriadjy)]<-sample.int(2,n*(n-1)/2,replace=TRUE)-1 # half filled

redocols<-which(colSums(lowertriadjy)>K)

for(j in redocols){
  currentparents<-which(lowertriadjy[,j]>0)
  lowertriadjy[sample(currentparents,length(currentparents)-K),j]<-0
}

# DAG data will be generated from

incidence<-adjpermtoDAG(n,nodelabels,lowertriadjy)

# permutation and partition data generated from

realpartypermy<-DAGtopartition(n,incidence) # turn the new DAG into a partition and permutation
realpermy<-realpartypermy$permy 
realparty<-realpartypermy$party

# generate some simulated Data

xarray<-matrix(0,n,N)

for(i in c(1:n)){
  currentnode<-realpermy[n-i+1]
  xtemp<- rnorm(N,mean=0,sd=sqrt(0.2))
  parents<-which(incidence[,currentnode]>0)
  numpars<-length(parents)
  if(numpars>0){
    for(j in parents){
	xtemp<-xtemp+2*xarray[j,]
    }
  }
  #xarray[currentnode,]<-(xtemp-mean(xtemp))/sd(xtemp) # rescale at each step
  xarray[currentnode,]<-xtemp # or do not rescale at each step?
}

Data <- xarray





