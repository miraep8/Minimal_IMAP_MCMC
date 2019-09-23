
# scoring parameters

am<-1
aw<-n+am+1
mu0<-numeric(n)
T0scale <- am*(aw-n-1)/(am+1) # This follows from equations (19) and (20) of [GH2002]
T0<-diag(T0scale,n,n)
TN <- T0 + (N-1)* cov(t(Data)) + ((am*N)/(am+N))* (mu0 - rowMeans(Data))%*%t(mu0 - rowMeans(Data))

awpN<-aw+N
constscorefact<- -(N/2)*log(pi) + (1/2)*log(am/(am+N))

scoreconstvec<-numeric(n)
for (j in 1:n){# j represents the number of parents plus 1
awp<-aw-n+j
scoreconstvec[j]<-constscorefact - lgamma(awp/2) + lgamma((awp+N)/2) + ((awp+j-1)/2)*log(T0scale)
}







