# This version calculates the BGe score using
# more numerically stable methods for when 
# the elements of TN become too disparate

# It is marginally slower than the previous version
# but the difference will be negligible for large graphs

# The determinant of a 2 by 2 matrix

dettwobytwo <- function(D){
  D[1,1]*D[2,2]-D[1,2]*D[2,1]
}

# The determinant of a 3 by 3 matrix

detthreebythree <- function(D){
	D[1,1]*(D[2,2]*D[3,3]-D[2,3]*D[3,2])-D[1,2]*(D[2,1]*D[3,3]-D[2,3]*D[3,1])+D[1,3]*(D[2,1]*D[2,3]-D[2,2]*D[3,1])
}

# The log of the BGe score, but simplified as much as possible
# see arXiv:1402.6863 

DAGcorescore<-function(j,parentnodes,n){

  lp<-length(parentnodes) #number of parents
  awpNd2<-(awpN-n+lp+1)/2
  
  A<-TN[j,j]
  
  switch(as.character(lp),
         "0"={# just a single term if no parents
           corescore <- scoreconstvec[lp+1] -awpNd2*log(A)
         },
         
         "1"={# no need for matrices
           D<-TN[parentnodes,parentnodes]
           logdetD<-log(D)
           B<-TN[j,parentnodes]
           logdetpart2<-log(A-B^2/D)
           corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
         },
         
         "2"={# can do matrix determinant and inverse explicitly
           # but this is numerically unstable for large matrices!
           # so we use the same approach as for 3 parents
           D<-TN[parentnodes,parentnodes]
           detD<-dettwobytwo(D)
           logdetD<-log(detD)
           B<-TN[j,parentnodes]
           #logdetpart2<-log(A-(D[2,2]*B[1]^2+D[1,1]*B[2]^2-2*D[1,2]*B[1]*B[2])/detD) #also using symmetry of D
           logdetpart2<-log(dettwobytwo(D-(B)%*%t(B)/A))+log(A)-logdetD
           corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
         },
         
         #"3"={# can still do matrix determinants efficiently so we use other approach
           # but the explicit formula is numerically unstable!
        #   D<-TN[parentnodes,parentnodes]
        #   detD<-detthreebythree(D)
        #   logdetD<-log(detD)
        #   B<-TN[j,parentnodes]
        #   logdetpart2<-log(detthreebythree(D-(B)%*%t(B)/A))+log(A)-logdetD
        #   corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
        # },
         
  {# otherwise we use cholesky decomposition to perform both
    D<-as.matrix(TN[parentnodes,parentnodes])
    choltemp<-chol(D)
    logdetD<-2*log(prod(choltemp[(lp+1)*c(0:(lp-1))+1]))
    B<-TN[j,parentnodes]
    logdetpart2<-log(A-sum(backsolve(choltemp,B,transpose=TRUE)^2))
    corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
  })
  
}
