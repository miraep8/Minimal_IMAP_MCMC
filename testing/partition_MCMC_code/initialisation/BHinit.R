
load("./initialisation/BostonHdata.Rdata")

# Size of the problem

Data<-BostonHdata

n <- nrow(Data) # number of nodes
N <- ncol(Data) # number of observations

# possible best fitting DAGs 

incidence<-matrix(0,n,n)
incidence[9,1]<-1
incidence[13,1]<-1
incidence[5,3]<-1
incidence[6,3]<-1
incidence[8,3]<-1
incidence[5,4]<-1
incidence[14,4]<-1
incidence[2,6]<-1
incidence[5,7]<-1
incidence[8,7]<-1
incidence[13,7]<-1
incidence[3,9]<-1
incidence[5,9]<-1
incidence[2,10]<-1
incidence[3,10]<-1
incidence[9,10]<-1
incidence[2,11]<-1
incidence[8,11]<-1
incidence[10,11]<-1
incidence[9,12]<-1
incidence[11,12]<-1
incidence[14,12]<-1
incidence[5,13]<-1
incidence[8,13]<-1
incidence[14,13]<-1
incidence[5,14]<-1
incidence[6,14]<-1
incidence[11,14]<-1

# both edges below could be reversed but this doesn't change the score

incidence[8,5]<-1
incidence[2,8]<-1

# possible best fitting permutations and partitions

realpermy<-c(1,7,4,12,13,14,11,10,9,3,5,6,8,2)
realparty<-c(2,3,1,1,1,1,1,1,2,1)

#realpermy<-c(1,7,4,12,13,14,11,10,9,3,6,2,8,5) # this corresponds to the reversed edges
#realparty<-c(2,3,1,1,1,1,1,1,1,1,1) # ditto
