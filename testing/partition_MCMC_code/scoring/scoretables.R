# This function produces a matrix with all the possible parents of a given node
# up to the maximum number of parents

listpossibleparents<-function(maxparents,elements){

  listy<-vector("list",length(elements))

  for (i in elements){
    remainingelements<-elements[-i]

    matrixofparents<-rep(NA,maxparents)
    for (r in 1:maxparents){
      possparents<-combinations(length(remainingelements),r,remainingelements)
      if(r<maxparents){
        for (j in 1:(maxparents-r)){
	  possparents <- cbind(possparents, NA)
        }
      }
    matrixofparents<-rbind(matrixofparents,possparents,deparse.level=0)
    }
  listy[[i]] <- matrixofparents
  }  

return(listy)
}

# This function scores all the possible parents

scorepossibleparents<-function(parenttable,n){

  listy<-vector("list",n)

  for (j in 1:n){
    scoretemp<-TableDAGscore(parenttable[[j]], j, n)
    listy[[j]] <- as.matrix(scoretemp)
  }  

return(listy)

}

# A compare function that treats NA as an element

compareNA <- function(v1,v2) {
    # This function returns TRUE wherever elements are the same, including NA's,
    # and false everywhere else.
    same <- (v1 == v2)  |  (is.na(v1) & is.na(v2))
    same[is.na(same)] <- FALSE
    return(same)
}

# This function scores certain nodes of an incidence matrix given the parent and score tables

DAGscorefromtable <- function(incidence, n, rescorenodes,parenttable,scoretable){
	P_local <- numeric(n)
	maxparents<-length(parenttable[[1]][1,])   
	for (j in rescorenodes)  {
		parentnodes <- which(incidence[,j]==1)
		parentline<-c(parentnodes,rep(NA,maxparents-length(parentnodes)))
		chosenline<-which(apply(parenttable[[j]],1,function(x) all(compareNA(x,parentline))))
		P_local[j]<-scoretable[[j]][chosenline]
	}
	return(P_local)
}
