
orderMCMC<-function(n,startorder,iterations,stepsave,parenttable,scoretable,moveprobs){

  currentpermy<-startorder #starting order represented as a permutation
  currentorderscores<-orderscore(n,c(1:n),parenttable,scoretable,currentpermy) #starting score
  currenttotallogscore<-sum(currentorderscores$totscores) #log total score of all DAGs in the starting order

  currentDAG<-samplescore(n,currentorderscores) #score of a single DAG sampled from the starting order

  L1 <- list() # stores the adjacency matrix of a DAG sampled from the orders
  L2 <- list() # stores its log BGe score
  L3 <- list() # stores the log BGe score of the entire order
  L4 <- list() # stores the orders as permutations

  zlimit<- floor(iterations/stepsave) + 1 # number of outer iterations
  length(L1) <- zlimit
  length(L2) <- zlimit
  length(L3) <- zlimit
  length(L4) <- zlimit

  L1[[1]]<-currentDAG$incidence #starting DAG adjacency matrix
  L2[[1]]<-currentDAG$logscore #starting DAG score
  L3[[1]]<-currenttotallogscore #starting order score
  L4[[1]]<-currentpermy #starting order

  for (z in 2:zlimit){ #the MCMC chain loop with 'iteration' steps is in two parts
    for (count in 1:stepsave){ #since we only save the results to the lists each 'stepsave'

      chosenmove<-sample.int(3,1,prob=moveprobs)
      if(chosenmove<3){	# if it is 3 then we stay still
	
	proposedpermy<-currentpermy #sample a new order by swapping two elements
	switch(as.character(chosenmove),
	"1"={ # swap any two elements at random
	  sampledelements<-sample.int(n,2,replace=FALSE) #chosen at random
	},
	"2"={ # swap any adjacent elements
	  k<-sample.int(n-1,1) #chose the smallest at random
	  sampledelements<-c(k,k+1)
	},
	{# if neither is chosen, we have a problem
	  print('The move sampling has failed!')
	})
	proposedpermy[sampledelements]<-currentpermy[rev(sampledelements)] #proposed new order

	rescorenodes<-proposedpermy[min(sampledelements):max(sampledelements)] #we only need to rescore these nodes between the swapped elements to speed up the calculation

	proposedorderrescored<-orderscore(n,rescorenodes,parenttable,scoretable,proposedpermy) #their scores
	proposedtotallogscore<-currenttotallogscore-sum(currentorderscores$totscores[rescorenodes])+sum(proposedorderrescored$totscores[rescorenodes]) #and the new log total score by updating only the necessary nodes

	scoreratio<-exp(proposedtotallogscore-currenttotallogscore) #acceptance probability

	if(runif(1)<scoreratio){ #Move accepted then set the current order and scores to the proposal
		currentpermy<-proposedpermy
		currentorderscores$allowedrows[rescorenodes]<-proposedorderrescored$allowedrows[rescorenodes]
		currentorderscores$allscores[rescorenodes]<-proposedorderrescored$allscores[rescorenodes]
		currentorderscores$totscores[rescorenodes]<-proposedorderrescored$totscores[rescorenodes]
		currenttotallogscore<-proposedtotallogscore
	}
      }
    }
    currentDAG<-samplescore(n,currentorderscores)
    L1[[z]]<-currentDAG$incidence #store adjacency matrix of a sampled DAG each 'stepsave'
    L2[[z]]<-currentDAG$logscore #and log score of a sampled DAG
    L3[[z]]<-currenttotallogscore #and the current order score
    L4[[z]]<-currentpermy #and store current order
  }
  return(list(L1,L2,L3,L4))
}

