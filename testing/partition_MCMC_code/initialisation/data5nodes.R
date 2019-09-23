set.seed(333)

N<-100 # Number of observations

### Example for the talk

# the commented out lines below would rescale the data at each step

x1 <- rnorm(N, mean=0, sd=sqrt(0.2))
#x1 <- (x1 - mean(x1))/sd(x1)

x3 <- rnorm(N, mean=0, sd=sqrt(0.2))
#x3 <- (x3 - mean(x3))/sd(x3)

x5 <- rnorm(N, mean=0, sd=sqrt(0.2))
#x5 <- (x5 - mean(x5))/sd(x5)

x4 <- 2*x3 + 2*x5 +  rnorm(N, mean=0, sd=sqrt(0.2))
#x4 <- (x4 - mean(x4))/sd(x4)

x2 <- 2*x1 + 2*x3 + 2*x4 + 2*x5 + rnorm(N, mean=0, sd=sqrt(0.2))
#x2 <- (x2 - mean(x2))/sd(x2)

Data <- rbind(x1, x2, x3, x4, x5)

# Size of the problem

n <- nrow(Data) # number of nodes
N <- ncol(Data) # number of observations

# DAG data generated from

incidence<-matrix(0,n,n)
incidence[5,4]<-1
incidence[5,2]<-1
incidence[4,2]<-1
incidence[3,4]<-1
incidence[3,2]<-1
incidence[1,2]<-1

# possible orders data generated from

realorderpermys<-list()
realorderpermys[[1]]<-c(2,4,1,3,5)
realorderpermys[[2]]<-c(2,4,1,5,3)
realorderpermys[[3]]<-c(2,4,3,1,5)
realorderpermys[[4]]<-c(2,4,3,5,1)
realorderpermys[[5]]<-c(2,4,5,1,3)
realorderpermys[[6]]<-c(2,4,5,3,1)
realorderpermys[[7]]<-c(2,1,4,3,5)
realorderpermys[[8]]<-c(2,1,4,5,3)

# permutation and partition data generated from

realpermy<-c(2,4,1,3,5)
realparty<-c(1,1,3)

