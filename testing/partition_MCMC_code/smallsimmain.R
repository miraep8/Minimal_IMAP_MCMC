# Choose the MCMC scheme

MCMCtype<-1 # 1 means standard structure, 2 with new edge reversal
# 3 means order MCMC, 4 means partition MCMC
# 5 means partition MCMC with new edge reversal

seedset<-1 # set a seed?
seednumbers<-100+c(1:2) # 101 and 102 are the values in the plots
seednumber<-seednumbers[1]

# Choose whether to save the plots to pdf
saveoutput<-1 # 1 means save

# load the necessary functions

source('./edgerevandstructure/structurefns.R')
source('./edgerevandstructure/structureMCMC.R')
source('./edgerevandstructure/newedgerevfns.R')
source('./edgerevandstructure/newedgerevmove.R')
source('./orderandpartition/orderMCMC.R')
source('./orderandpartition/orderfns.R')
source('./orderandpartition/partitionMCMC.R')
source('./orderandpartition/partitionmoves.R')
source('./orderandpartition/partitionfns.R')
source('./orderandpartition/samplefns.R')
source('./scoring/combinations.R')
source('./scoring/scorefns.R')
source('./scoring/scoretables.R')

# Here we use the BGe score

source('./scoring/bgescorefast.R')

# Generate data and load score parameters

source('./initialisation/data5nodes.R')
source('./initialisation/scoreparas.R')

# Score of the DAG data generated from

realDAGlogscore<-sum(DAGnodescore(incidence, n, c(1:n)))

# Choose maximum number of parents

maxparents<-4 # Maximum number of parents to allow

starttime<-proc.time() # for timing the problem

# Fill up a matrix with possible parents

parenttable<-listpossibleparents(maxparents,c(1:n))
tablelength<-nrow(parenttable[[1]]) # size of the table

# Now need to score them!

scoretable<-scorepossibleparents(parenttable, n) 

endtime<-proc.time()
endtime<-endtime-starttime
print('Time to initialise the score table')
print(endtime)

# Start the main part

switch(as.character(MCMCtype),
  "1"={ # standard structure MCMC
	iterations<-50000 #number of iterations in the chain
	moveprobs<-c(1) # having length 1 disallows the new edge reversal move
	if(!(length(moveprobs)==1)){print('Vector of move probabilities has the wrong length!')}
  },
  "2"={ # with new edge reversal
	iterations<-40000 #number of iterations in the chain
# Choose the probability of the different moves
# 1 is standard structure MCMC (including the possibility to stay still [officially needed for convergence])
# 2 is new edge reversal move
	moveprobs<-c(0.93,0.07)
	moveprobs<-moveprobs/sum(moveprobs) # normalisation
	if(!(length(moveprobs)==2)){print('Vector of move probabilities has the wrong length!')}
  },
  "3"={ # order MCMC
	iterations<-20000 #number of iterations in the chain
# Choose the probability of the different moves
# 1 is swap any two elements
# 2 is to only swap adjacent elements
# 3 is to stay still (officially needed for convergence)
	prob1<-99
	if(n>3){ prob1<-round(6*99*n/(n^2+10*n-24)) }
	prob1<-prob1/100
	moveprobs<-c(prob1,0.99-prob1,0.01)
	moveprobs<-moveprobs/sum(moveprobs) # normalisation
	if(!(length(moveprobs)==3)){print('Vector of move probabilities has the wrong length!')}
# calculate score of orders compatible with the DAG that the data is generated from
	realorderlogscores<-rep(0,length(realorderpermys))
	for (j in 1:length(realorderpermys)){
	realorderscores<-orderscore(n,c(1:n),parenttable,scoretable,realorderpermys[[j]])
	realorderlogscores[j]<-sum(realorderscores$totscores) #log total score of all DAGs in the order
	}
  },
  "4"={ # partition MCMC
	iterations<-10000 #number of iterations in the chain
# Choose the probability of the different moves
# 1 is swap two nodes from different partition elements
# 2 is to only swap nodes from adjacent elements
# 3 is to split or join partition elements
# 4 is to move a single node
# 5 is to stay still (officially needed for convergence)
	prob1start<-40/100
	prob1<-prob1start*100
	if(n>3){ prob1<-round(6*prob1*n/(n^2+10*n-24)) }
	prob1<-prob1/100
	prob2start<-99/100-prob1start
	prob2<-prob2start*100
	if(n>3){ prob2<-round(6*prob2*n/(n^2+10*n-24)) }
	prob2<-prob2/100
	moveprobs<-c(prob1,prob1start-prob1,prob2start-prob2,prob2,0.01)
	moveprobs<-moveprobs/sum(moveprobs) # normalisation
	if(!(length(moveprobs)==5)){print('Vector of move probabilities has the wrong length!')}
# calculate score of partition and permutation that the data is generated from
	realposy<-parttolist(n,realparty)
	realpartitionscores<-partitionscore(n,c(1:n),parenttable,scoretable,realpermy,realparty,realposy)
	realpartitionlogscore<-sum(realpartitionscores$totscores) #log total score of all DAGs in the partition
  },
  "5"={ # partition MCMC with edge reversal
	iterations<-9000 #number of iterations in the chain
# Choose the probability of the different moves
# 1 is swap two nodes from different partition elements
# 2 is to only swap nodes from adjacent elements
# 3 is to split or join partition elements
# 4 is to move a single node
# 5 is to stay still (officially needed for convergence)
# 6 is the new edge reversal
	prob1start<-37/100
	prob1<-prob1start*100
	if(n>3){ prob1<-round(6*prob1*n/(n^2+10*n-24)) }
	prob1<-prob1/100
	prob2start<-92/100-prob1start
	prob2<-prob2start*100
	if(n>3){ prob2<-round(6*prob2*n/(n^2+10*n-24)) }
	prob2<-prob2/100
	moveprobs<-c(prob1,prob1start-prob1,prob2start-prob2,prob2,0.01,0.07)
	moveprobs<-moveprobs/sum(moveprobs) # normalisation
	if(!(length(moveprobs)==6)){print('Vector of move probabilities has the wrong length!')}
# calculate score of partition and permutation that the data is generated from
	realposy<-parttolist(n,realparty)
	realpartitionscores<-partitionscore(n,c(1:n),parenttable,scoretable,realpermy,realparty,realposy)
	realpartitionlogscore<-sum(realpartitionscores$totscores) #log total score of all DAGs in the partition
  },
  {# if none is chosen, we have a problem
    print('Not implemented')
  })

# Plotting parameters

par(mfrow=c(2,1))

par(mar=c(2.5,5.75,0.5,0.75))
par(mgp=c(3.5,1.25,0))
par(cex.axis=1.25)
par(cex.lab=1.5)

stepsave<-iterations/1000 #how often to save the result

if(seedset>0){
  set.seed(seednumber) # choose one?
}

if(MCMCtype<3){

  startDAG<-matrix(numeric(n*n),nrow=n) # starting DAG is empty say

starttime<-proc.time() # for timing the problem

  revallowed<-1 # allow standard edge reversals

  example<-structureMCMC(n,startDAG,iterations,stepsave,maxparents,parenttable,scoretable,revallowed,moveprobs)

endtime<-proc.time()
endtime<-endtime-starttime
print('Time to run the')
print(iterations)
print('iterations')
print(endtime)

  revallowed<-0 # don't allow edge reversals

  if(seedset>0){
    set.seed(seednumber)
  }

  startDAG<-matrix(numeric(n*n),nrow=n) # starting DAG is empty say

starttime<-proc.time() # for timing the problem

  example2<-structureMCMC(n,startDAG,iterations,stepsave,maxparents,parenttable,scoretable,revallowed,moveprobs)

endtime<-proc.time()
endtime<-endtime-starttime
print('Time to run the')
print(iterations)
print('iterations')
print(endtime)

# Plot the results

  nparts<-length(example[[2]])
  maxDAGscore<-max(unlist(example[[2]]))

  plot(1:nparts,example[[2]],type="l", ylab="DAG score", xlab="",#xlab="iteration step", 
main="", col="blue",ylim=c(maxDAGscore-6.8,maxDAGscore+0.2))
  abline(h=maxDAGscore,col='springgreen4', lwd=3)
  abline(h=realDAGlogscore,col='red',lty=3,lwd=3)
  lines(1:nparts,example[[2]],type="l", col="blue")

  maxDAGscore2<-max(unlist(example2[[2]]))

  plot(1:nparts,example2[[2]],type="l", ylab="DAG score", xlab="",#xlab="iteration step", 
main="", col="blue",ylim=c(maxDAGscore2-6.8,maxDAGscore2+0.2))
  abline(h=maxDAGscore2,col='springgreen4', lwd=3)
  abline(h=realDAGlogscore,col='red',lty=3, lwd=3)
  lines(1:nparts,example2[[2]],type="l", col="blue")

} else if (MCMCtype==3){

  startorder<-c(1:n) # starting order

starttime<-proc.time() # for timing the problem

  example<-orderMCMC(n,startorder,iterations,stepsave,parenttable,scoretable,moveprobs)

endtime<-proc.time()
endtime<-endtime-starttime
print('Time to run the')
print(iterations)
print('iterations')
print(endtime)

# Plot the results

  nparts<-length(example[[2]])
  maxorderscore<-max(unlist(example[[3]]))
  maxDAGscore<-max(unlist(example[[2]]))

  plot(1:nparts,example[[3]],type="l", ylab="Order score", xlab="",#xlab="iteration step", 
main="", col="blue",ylim=c(maxorderscore-6.8,maxorderscore+0.2))
  abline(h=maxorderscore,col='springgreen4', lwd=3)
  abline(h=realorderlogscores,col='red',lty=3,lwd=3)
  lines(1:nparts,example[[3]],type="l", col="blue")

  plot(1:nparts,example[[2]],type="l", ylab="DAG score", xlab="",#xlab="iteration step", 
main="", col="blue",ylim=c(maxDAGscore-6.8,maxDAGscore+0.2))
  abline(h=maxDAGscore,col='springgreen4', lwd=3)
  abline(h=realDAGlogscore,col='red',lty=3,lwd=3)
  lines(1:nparts,example[[2]],type="l", col="blue")
} else if (MCMCtype>3){

  startpermutation<-c(1:n) # pick a starting permutation
  startpartition<-c(n) # and a starting partition - c(n) gives the empty DAG

starttime<-proc.time() # for timing the problem

  example<-partitionMCMC(n,startpermutation,startpartition,iterations,stepsave,parenttable,scoretable,moveprobs)

endtime<-proc.time()
endtime<-endtime-starttime
print('Time to run the')
print(iterations)
print('iterations')
print(endtime)

# Plot the results

  nparts<-length(example[[2]])
  maxpartitionscore<-max(unlist(example[[3]]))
  maxDAGscore<-max(unlist(example[[2]]))

  plot(1:nparts,example[[3]],type="l", ylab="Partition score", xlab="",#xlab="iteration step", 
main="", col="blue",ylim=c(maxpartitionscore-6.8,maxpartitionscore+0.2))
  abline(h=maxpartitionscore,col='springgreen4', lwd=3)
  abline(h=realpartitionlogscore,col='red',lty=3,lwd=3)
  lines(1:nparts,example[[3]],type="l", col="blue")

  plot(1:nparts,example[[2]],type="l", ylab="DAG score", xlab="",#xlab="iteration step", 
main="", col="blue",ylim=c(maxDAGscore-6.8,maxDAGscore+0.2))
  abline(h=maxDAGscore,col='springgreen4', lwd=3)
  abline(h=realDAGlogscore,col='red',lty=3,lwd=3)
  lines(1:nparts,example[[2]],type="l", col="blue")
}

if((saveoutput==1)&&(seedset==1)){
  switch(as.character(MCMCtype),
    "1"={ # standard structure MCMC
	dev.copy(pdf,paste("./smallsimgraphs/structure5nodes",N,"its",iterations/1000,"seed",seednumber,".pdf",sep=""), width=7.5, height=7.5, onefile=F, pointsize=10,  paper="special")
	dev.off()
    },
    "2"={ # with new edge reversal
	dev.copy(pdf,paste("./smallsimgraphs/edgerev5nodes",N,"one",100*moveprobs[1],"its",iterations/1000,"seed",seednumber,".pdf",sep=""), width=7.5, height=7.5, onefile=F, pointsize=10,  paper="special")
	dev.off()
    },
    "3"={ # order MCMC
	dev.copy(pdf,paste("./smallsimgraphs/order5nodes",N,"one",100*moveprobs[1],"two",100*moveprobs[2],"its",iterations/1000,"seed",seednumber,".pdf",sep=""), width=7.5, height=7.5, onefile=F, pointsize=10,  paper="special")
	dev.off()
    },
    "4"={ # Partition MCMC
	dev.copy(pdf,paste("./smallsimgraphs/partition5nodes",N,"one",100*moveprobs[1],"two",100*moveprobs[2],"three",100*moveprobs[3],"its",iterations/1000,"seed",seednumber,".pdf",sep=""), width=7.5, height=7.5, onefile=F, pointsize=10,  paper="special")
	dev.off()
    },
    "5"={ # Partition MCMC with edge reversal
	dev.copy(pdf,paste("./smallsimgraphs/partitionedgerev5nodes",N,"one",100*moveprobs[1],"two",100*moveprobs[2],"three",100*moveprobs[3],"six",100*moveprobs[6],"its",iterations/1000,"seed",seednumber,".pdf",sep=""), width=7.5, height=7.5, onefile=F, pointsize=10,  paper="special")
	dev.off()
    },
    {# if none is chosen, we have a problem
    print('Not implemented')
    })
}
