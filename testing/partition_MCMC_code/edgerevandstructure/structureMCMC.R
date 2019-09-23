# this code is taken from the Dortmund course programmed by Miriam Lohr
# I have only changed the scoring function to be correct
# replaced sample by 'propersample'
# and corrected the fan.in for edge reversal!

# The possibility of the new edge reversal move is also built in
# It is activated as long as moveprobs has two components
# The second being the probability of choosing the new edge reversal

# Standard edge reversal can also be toggled with 'revallowed'

structureMCMC <- function(n,incidence,iterations,stepsave,fan.in,parenttable,scoretable,revallowed,moveprobs){

  newedgerevallowed<-1 # allow the new edge reversal moves
  if(length(moveprobs)==1){
    newedgerevallowed<-0 # or don't allow them
    chosenmove<-1 
  }

  currentDAGlogscores<-DAGnodescore(incidence, n, c(1:n)) # score directly
  #currentDAGlogscores<-DAGscorefromtable(incidence,n,c(1:n),parenttable,scoretable) # or from the table of scores
  currentDAGlogscore<-sum(currentDAGlogscores) # score of incidence matrix

  L1 <- list() #stores the adjecency matrices
  L2 <- list() # stores the log BGe score of the DAGs

  zlimit<- floor(iterations/stepsave) + 1 # number of outer iterations
  length(L1) <- zlimit
  length(L2) <- zlimit

  L1[[1]]<-incidence #starting adjacency matrix
  L2[[1]]<-currentDAGlogscore #starting DAG score

# first ancestor matrix
  ancest1 <- ancestor(incidence)

####### ... the number of neighbour graphs/proposal probability for the FIRST graph
### 1.) number of neighbour graphs obtained by edge deletions
  num_deletion <- sum(incidence)

  emptymatrix<-matrix(numeric(n*n),nrow=n)
  fullmatrix<-matrix(rep(1,n*n),nrow=n)
  Ematrix<-diag(1,n,n)

### 2.) number of neighbour graphs obtained by edge additions    1- E(i,j) - I(i,j) - A(i,j)
  inter_add <- which(fullmatrix - Ematrix - incidence - ancest1 >0)
  add <- emptymatrix
  add[inter_add] <- 1
  add[,which(colSums(incidence)>fan.in-1)] <- 0
  num_addition <- sum(add)

### 3.) number of neighbour graphs obtained by edge reversals    I - (I^t * A)^t
  inter_rev <- which(incidence - t(t(incidence)%*% ancest1)==1)
  re <- emptymatrix
  re[inter_rev] <- 1
  re[which(colSums(incidence)>fan.in-1),] <- 0 # this has to be this way around
  num_reversal <- sum(re)

##### total number of neighbour graphs:
  currentnbhoodnorevs <- sum(num_deletion,num_addition)+1
  currentnbhood <- currentnbhoodnorevs+num_reversal

  for (z in 2:zlimit){
    for (count in 1:stepsave){

      if(newedgerevallowed>0){ # if we allow the new edge reversal move then sample the move type	
        chosenmove<-sample.int(2,1,prob=moveprobs) # sample what type of move
      }
      switch(as.character(chosenmove),
      "1"={ # standard structure MCMC

### sample one of the three single edge operations
	if(revallowed==1){
	  operation<-sample.int(4,1,prob=c(num_reversal,num_deletion,num_addition,1)) # sample the type of move including staying still
	} else {
	  operation<-sample.int(3,1,prob=c(num_deletion,num_addition,1))+1 # sample the type of move including staying still
	}

# 1 is edge reversal, 2 is deletion and 3 is additon. 4 represents choosing the current DAG

  if(operation<4){ # if we don't stay still

#### shifting of the incidence matrix
	incidence_new <- incidence

# creating a matrix with dimensions of the incidence matrix and all entries zero except for the entry of the chosen edge
	help_matrix <- emptymatrix

	if (operation==2){              # if edge deletion was sampled
	  new_edge <- propersample(which(incidence>0)) # sample one of the existing edges
	  incidence_new[new_edge] <- 0            # and delete it
	  help_matrix[new_edge] <- 1
	}
	if (operation==1){    # if edge reversal was sampled
	  new_edge <- propersample(which(re==1))      # sample one of the existing edges where a reversal leads to a valid graph
	  incidence_new[new_edge] <- 0             # delete it
  	  help_matrix[new_edge] <- 1               # an only a "1" at the entry of the new (reversed) edge
	  incidence_new <- incidence_new + t(help_matrix) # sum the deleted matrix and the "help-matrix"
	}
	if (operation==3){     # if edge addition was sampled
	  new_edge <- propersample(which(add==1)) # sample one of the existing edges where a addition leads to a valid graph
	  incidence_new[new_edge] <- 1             # and add it
	  help_matrix[new_edge] <- 1
	}

### Updating the ancestor matrix

# numbers of the nodes that belong to the shifted egde
	parent <- which(rowSums(help_matrix)==1)
	child <- which(colSums(help_matrix)==1)

### updating the ancestor matrix (after edge reversal)
## edge deletion
	ancestor_new <- ancest1
	if (operation==1){                                                                          
	  ancestor_new[c(child,which(ancest1[,child]==1)),] <- 0   # delete all ancestors of the child and its descendants                                           
	  top_name <- des_top_order(incidence_new, ancest1, child)
	  for (d in top_name){
	    for(g in which(incidence_new[,d]==1)) {
		ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
	    }
	  }

	  anc_parent <- which(ancestor_new[child,]==1)          # ancestors of the new parent
	  des_child <- which(ancestor_new[,parent]==1)          # descendants of the child
	  ancestor_new[c(parent,des_child),c(child,anc_parent)] <- 1
	}

### updating the ancestor matrix (after edge deletion)
	if (operation==2){
	  ancestor_new[c(child,which(ancest1[,child]==1)),] <- 0   # delete all ancestors of the child and its descendants                                           #
	  top_name <- des_top_order(incidence_new, ancest1, child)
	  for (d in top_name){
	    for(g in which(incidence_new[,d]==1)) {
		ancestor_new[d,c(g,(which(ancestor_new[g,]==1)))] <- 1
	    }
	  }
	}

# updating the ancestor matrix (after edge addition)
	if (operation==3){
	  anc_parent <- which(ancest1[parent,]==1)             # ancestors of the new parent
	  des_child <- which(ancest1[,child]==1)               # descendants of the child
	  ancestor_new[c(child,des_child),c(parent,anc_parent)] <- 1
	}

####### ... the number of neighbour graphs/proposal probability for the proposed graph
### 1.) number of neighbour graphs obtained by edge deletions
	num_deletion_new <- sum(incidence_new)

### number of neighbour graphs obtained by edge additions    1- E(i,j) - I(i,j) - A(i,j)
	inter_add.new <- which(fullmatrix - Ematrix - incidence_new - ancestor_new >0)
	add.new <- emptymatrix
	add.new[inter_add.new] <- 1
	add.new[,which(colSums(incidence_new)>fan.in-1)] <- 0
	num_addition_new <- sum(add.new)

### number of neighbour graphs obtained by edge reversals    I - (I^t * A)^t
	inter_rev.new<- which(incidence_new - t(t(incidence_new)%*% ancestor_new)==1)
	re.new <- emptymatrix
	re.new[inter_rev.new] <- 1
	re.new[which(colSums(incidence_new)>fan.in-1),] <- 0 # this has to be this way around!
	num_reversal_new <- sum(re.new)

##### total number of neighbour graphs:
	proposednbhoodnorevs<-sum(num_deletion_new, num_addition_new) + 1
	proposednbhood <- proposednbhoodnorevs + num_reversal_new

	rescorenodes <- child

	if (operation==1){                      # if single edge operation was an edge reversal
	  rescorenodes<-c(child,parent)
	}

  proposedDAGrescored<-DAGnodescore(incidence_new, n, rescorenodes) # rescore relevant nodes directly
	#proposedDAGrescored<-DAGscorefromtable(incidence_new, n, rescorenodes,parenttable,scoretable) # or from the score table
	proposedDAGlogscore<-currentDAGlogscore-sum(currentDAGlogscores[rescorenodes])+sum(proposedDAGrescored[rescorenodes]) #and the new log total score by updating only the necessary nodes

	if(revallowed==1){
	  scoreratio<-exp(proposedDAGlogscore-currentDAGlogscore)*(currentnbhood/proposednbhood) #acceptance probability
	} else{
	  scoreratio<-exp(proposedDAGlogscore-currentDAGlogscore)*(currentnbhoodnorevs/proposednbhoodnorevs) #acceptance probability
	}

	if(runif(1)<scoreratio){ #Move accepted then set the current order and scores to the proposal
	  incidence <- incidence_new
	  currentDAGlogscores[rescorenodes]<-proposedDAGrescored[rescorenodes]
	  currentDAGlogscore<-proposedDAGlogscore
	  ancest1 <- ancestor_new
	  currentnbhood<-proposednbhood
	  currentnbhoodnorevs<-proposednbhoodnorevs	
	  num_deletion <- num_deletion_new
	  num_addition <- num_addition_new
	  num_reversal <- num_reversal_new
	  add <- add.new
	  re <- re.new
	}
  }  # end of staying still loop
      },
      "2"={ # new edge reversal move
	testy<-newedgereversalmove(n,incidence,parenttable,scoretable)
	if(is.list(testy)){ # if move was accepted
	  incidence<-testy$incidence # update the matrix
    rescorenodes<-testy$rescore # the child and parent which were swapped
	  currentDAGlogscores[rescorenodes]<-DAGnodescore(incidence, n, rescorenodes)[rescorenodes]
	  currentDAGlogscore<-sum(currentDAGlogscores) # score of incidence matrix 

# The new ancestor matrix
    ancest1<-testy$ancestor

####### ... the number of neighbour graphs/proposal probability for the FIRST graph
### 1.) number of neighbour graphs obtained by edge deletions
	  num_deletion <- sum(incidence)

### 2.) number of neighbour graphs obtained by edge additions    1- E(i,j) - I(i,j) - A(i,j)
	  inter_add <- which(fullmatrix - Ematrix - incidence - ancest1 >0)
	  add <- emptymatrix
	  add[inter_add] <- 1
	  add[,which(colSums(incidence)>fan.in-1)] <- 0
	  num_addition <- sum(add)

### 3.) number of neighbour graphs obtained by edge reversals    I - (I^t * A)^t
	  inter_rev <- which(incidence - t(t(incidence)%*% ancest1)==1)
	  re <- emptymatrix
	  re[inter_rev] <- 1
	  re[which(colSums(incidence)>fan.in-1),] <- 0 # and correct here too!
	  num_reversal <- sum(re)

##### total number of neighbour graphs:
	  currentnbhoodnorevs <- sum(num_deletion,num_addition) + 1
	  currentnbhood <- currentnbhoodnorevs + num_reversal
	}
      },
      {# if neither is chosen, we have a problem
	    print('The move sampling has failed!')
      })
    }
	
  L1[[z]] <- incidence
  L2[[z]] <- currentDAGlogscore 
  }
return(list(L1,L2))
}

################################################################################

