
# This function samples an element from a vector properly

propersample <- function(x){if(length(x)==1) x else sample(x,1)} 

# This function samples a single DAG score according to a list of possible log scores

samplescore<-function(n,scores){
  incidence<-matrix(numeric(n*n),nrow=n) # store the adjacency matrix
  sampledscore<-0
  for (i in 1:n){
    scorelength<-length(scores$allscores[[i]])
    k<-sample.int(scorelength,1,prob=exp(scores$allscores[[i]]-scores$totscores[i])) # sample according to scores
    parentrow<-parenttable[[i]][scores$allowedrows[[i]][k],] # the parent set
    parentset<-parentrow[which(parentrow>0)] # removing NAs
    incidence[parentset,i]<-1 # fill in elements of the adjacency matrix
    sampledscore<-sampledscore+scores$allscores[[i]][k] # and add the score
  }
  DAG<-list()
  DAG$incidence<-incidence
  DAG$logscore<-sampledscore
  return(DAG)
}
