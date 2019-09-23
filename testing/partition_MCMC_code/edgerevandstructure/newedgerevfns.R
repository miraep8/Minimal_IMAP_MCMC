
# This function gives the list of row labels of parentsets compatible with 
# child as a new parent and no banned nodes

possibleparentsedgerev<-function(n,parenttableentry,bannednodes,child){

  tablesize<-dim(parenttableentry) # just to remove some arguments

  allowedrows<-c(2:tablesize[1])
  for (j in 1:tablesize[2]){ # working columnwise allows R to speed up 
    bannedrows<-which(parenttableentry[allowedrows,j]%in%bannednodes)
    if(length(bannedrows)>0){
      allowedrows<-allowedrows[-bannedrows]
    }
  }
  notrequiredrows<-allowedrows
  for (j in 1:tablesize[2]){ # now we remove the allowable rows instead
    requiredrows<-which(parenttableentry[notrequiredrows,j]%in%child)
    if(length(requiredrows)>0){
	notrequiredrows<-notrequiredrows[-requiredrows]
    }
  }
  allowedrows<-setdiff(allowedrows,notrequiredrows) # and keep just the difference!

return(allowedrows)
}

# This function gives the list of row lables of parentsets 
# containing no banned nodes

possibleparentsedgerevnext<-function(n,parenttableentry,bannednodes){

  tablesize<-dim(parenttableentry) # just to remove some arguments

  allowedrows<-c(2:tablesize[1])
  for (j in 1:tablesize[2]){ # working columnwise allows R to speed up 
    bannedrows<-which(parenttableentry[allowedrows,j]%in%bannednodes)
    if(length(bannedrows)>0){
      allowedrows<-allowedrows[-bannedrows]
    }
  }
  allowedrows<-c(1,allowedrows) # this row is always allowed

return(allowedrows)
}

### calculation of a descendents matrix:
descendents <- function(incidence){
  incidence1 <- t(incidence)
  incidence2 <- t(incidence)
  k <- 1
  while (k < nrow(incidence)){
    incidence1 <- incidence1%*%t(incidence)
    incidence2 <- incidence2 + incidence1
    k <-k+1
  }
  incidence2[which(incidence2[,]>0)] <- 1
  return(t(incidence2))
}

