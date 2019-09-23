
# This function takes in the adjacency matrix of a DAG and performs a
# new edge reversal move of Grzegorczyk and Husmeier.

# If the move is accepted it returns the new DAG, ancestor matrix and which nodes were chosen

newedgereversalmove<-function(n,incidence,parenttable,scoretable){

  tobereturned<-NA
  edges<-which(incidence==1)

  Ndagger<-length(edges) # store the number of edges

  if(Ndagger>0){ # if there is at least one edge to reverse!
     samplededge<-propersample(edges)

# creating a matrix with dimensions of the incidence matrix and all entries zero except for the entry of the chosen edge
    helpmatrix <- matrix(0,n,n)
    helpmatrix[samplededge] <- 1

# labels of the node that belong to the selected egde
    edgeparent <- which(rowSums(helpmatrix)==1)
    edgechild <- which(colSums(helpmatrix)==1)

# First we need to orphan the selected edge nodes

    orphan<-incidence
    orphan[,edgechild]<-0 # we orphan the child
    orphan[,edgeparent]<-0 # and the parent

# Here we calculate the descendents in a naive way - could be performed more efficiently

    descendorph<-descendents(orphan)

# Then from the orphan we need to make two steps
# First we find new parents for the old parent node
# It cannot include its descendents and must include the edge's child

    edgeparentbannednodes<-which(descendorph[edgeparent,]>0) 

# The rows of the parent table which satisfy the conditions

    possibleparentrows<-possibleparentsedgerev(n,parenttable[[edgeparent]],edgeparentbannednodes,edgechild)

# and their scores

    possibleparentscores<-scoretable[[edgeparent]][possibleparentrows]
    scorelength<-length(possibleparentscores)

# to sample and sum exponentials properly
    maxscore<-max(possibleparentscores)
    expparentscores<-exp(possibleparentscores-maxscore)

# sample a row accordingly

    sampledelement<-sample.int(scorelength,1,prob=expparentscores)

# Store the partition function

    Zstarparent<-sum(expparentscores)
    Zstarparentlogscale<-maxscore

# and find the parents

    newparents<-parenttable[[edgeparent]][possibleparentrows[sampledelement],]
    newparents<-newparents[which(newparents>0)]#remove the NAs

# Fill up the new adjacency matrix

    incidence1<-orphan
    incidence1[newparents,edgeparent]<-1 # add the new parents to the old parent

# update the descendent matrix

    descend1<-descendorph
    descend1[newparents,edgeparent]<-1 # add the new edges to the descendent matrix
    oldparentdesc<-which(descend1[edgeparent,]>0) # the descendents of the old parent
    descend1[newparents,oldparentdesc]<-1 # and edges to the parents descendents
    if(length(newparents)>1){ # find the ancestors of the parents 
      newancestors<-which(rowSums(descend1[,newparents])>0) 
    } else {
      newancestors<-which(descend1[,newparents]>0) 
    }
    descend1[newancestors,c(edgeparent,oldparentdesc)]<-1 # and edges from the new parents' ancestors to the old parent and its descendents

#descend12<-descendents(incidence1) # This was to check the update works
#print(sum((descend1-descend12)^2))

# Now we need to sample new parents for the previous child
# none of the current descendents are permissible

    edgechildbannednodes<-which(descend1[edgechild,]>0)

# The rows of the parent table which satisfy the conditions

    possibleparentrows<-possibleparentsedgerevnext(n,parenttable[[edgechild]],edgechildbannednodes)

# and their scores

    possibleparentscores<-scoretable[[edgechild]][possibleparentrows]
    scorelength<-length(possibleparentscores)

# to sample and sum exponentials properly
    maxscore<-max(possibleparentscores)
    expparentscores<-exp(possibleparentscores-maxscore)

# sample a row accordingly

    sampledelement<-sample.int(scorelength,1,prob=expparentscores)

# Store the partition function

    Zpluschild<-sum(expparentscores)
    Zpluschildlogscale<-maxscore

# and find the parents

    newparents<-parenttable[[edgechild]][possibleparentrows[sampledelement],]
    newparents<-newparents[which(newparents>0)]#remove the NAs

    if(length(newparents)>0){
# Update the new adjacency matrix

      incidence1[newparents,edgechild]<-1 # add the new parents to the old child

# Update the descendent matrix too
      descend1[newparents,edgechild]<-1 # add the new edges to the descendent matrix
      oldchilddesc<-which(descend1[edgechild,]>0) # the descendents of the old child
      descend1[newparents,oldchilddesc]<-1 # and edges to the childs descendents

      if(length(newparents)>1){ # find the ancestors of the parents 
        newancestors<-which(rowSums(descend1[,newparents])>0) 
      } else {
        newancestors<-which(descend1[,newparents]>0) 
      }
      descend1[newancestors,c(edgechild,oldchilddesc)]<-1 # and edges from the new parents' ancestors to the old child and its descendents
    }

#descend12<-descendents(incidence1) # This was to check the update works
#print(sum((descend1-descend12)^2))

    Ndaggertilde<-length(which(incidence1==1)) # store the new number of edges

# Next we need to make the same steps from the orphan but keeping the child and parent relation as before!
# This is solely needed for the acceptance ratio.

# First we find different parents for the child node
# It cannot include its descendents and must include its previous edge parent

    edgechildbannednodes<-which(descendorph[edgechild,]>0) 

# The rows of the parent table which satisfy the conditions

    possibleparentrows<-possibleparentsedgerev(n,parenttable[[edgechild]],edgechildbannednodes,edgeparent)

# and their scores

    possibleparentscores<-scoretable[[edgechild]][possibleparentrows]
    scorelength<-length(possibleparentscores)

# to sample and sum exponentials properly
    maxscore<-max(possibleparentscores)
    expparentscores<-exp(possibleparentscores-maxscore)

# sample a row accordingly

    sampledelement<-sample.int(scorelength,1,prob=expparentscores)

# Store the partition function

    Zstarchild<-sum(expparentscores)
    Zstarchildlogscale<-maxscore

# and find the parents

    newparents<-parenttable[[edgechild]][possibleparentrows[sampledelement],]
    newparents<-newparents[which(newparents>0)]#remove the NAs

# Now we don't actually need to fill up the new adjacency matrix

#incidence2<-orphan
#incidence2[newparents,edgechild]<-1 # add the new parents to the old parent

# update the descendent matrix

    descend2<-descendorph
    descend2[newparents,edgechild]<-1 # add the new edges to the descendent matrix
    childdesc<-which(descend2[edgechild,]>0) # descendendents of the child
    descend2[newparents,childdesc]<-1 # and their descendents
    if(length(newparents)>1){ # find the ancestors of the parents 
      newancestors<-which(rowSums(descend2[,newparents])>0) 
    } else {
      newancestors<-which(descend2[,newparents]>0) 
    }
    descend2[newancestors,c(edgechild,childdesc)]<-1 # and edges from the new parents' ancestors to the old parent and its descendents

#descend22<-descendents(incidence2) # This was to check the update works
#print(sum((descend2-descend22)^2))

# Now we need to sample new parents for the previous parent
# none of the current descendents are permissible

    edgeparentbannednodes<-which(descend2[edgeparent,]>0)

# The rows of the parent table which satisfy the conditions

    possibleparentrows<-possibleparentsedgerevnext(n,parenttable[[edgeparent]],edgeparentbannednodes)

# and their scores

    possibleparentscores<-scoretable[[edgeparent]][possibleparentrows]

# to sum exponentials properly
    maxscore<-max(possibleparentscores)
    expparentscores<-exp(possibleparentscores-maxscore)

# Store the partition function

    Zplusparent<-sum(expparentscores)
    Zplusparentlogscale<-maxscore

# Now we can finally calculate the accepance ratio

    scoreratio<-(Ndagger/Ndaggertilde)*(Zstarparent/Zstarchild)*(Zpluschild/Zplusparent)*exp(Zstarparentlogscale-Zstarchildlogscale+Zpluschildlogscale-Zplusparentlogscale)

# MH move

    if(runif(1)<scoreratio){ # if move accepted return the new DAG
	tobereturned<-list("incidence"=incidence1,"ancestor"=t(descend1),"rescore"=c(edgechild,edgeparent))
    }
  } else{ # if no edges to be sampled we need to stay still
  #print('no edges to sample - stay still')
  }

return(tobereturned)

}


