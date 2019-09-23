# First choose the MCMC scheme

MCMCtype<-1 # 1 means standard structure, 2 with new edge reversal
# 3 means order MCMC, 4 means partition MCMC
# 5 means partition MCMC with new edge reversal

seedset<-1 # set a seed?
seednumber<-101 # an example used

# Choose whether to save the plots to png or pdf
saveoutput<-1 # 1 means save
filetypes<-c(".png",".pdf") # choose png or pdf for the figures
filetypey<-1 # we use png here  

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

# load a simple score proportional to the number of edges in the DAG

source('./scoring/numedgescore.R')

# where to save the data

dirname<-"./simpleexamplegraphs/"

# Choose size of DAGs to consider

n<-3

# Choose maximum number of parents

maxparents<-2 # Maximum number of parents to allow

starttime<-proc.time() # for timing the problem

# Fill up a matrix with possible parents

parenttable<-listpossibleparents(maxparents,c(1:n))
tablelength<-nrow(parenttable[[1]]) # size of the table

# Now need to score them!

scoretable<-scorepossibleparents(parenttable,n) 

endtime<-proc.time()
endtime<-endtime-starttime
print('Time to initialise the score table')
print(endtime)

switch(as.character(MCMCtype),
  "1"={ # standard structure MCMC
	iterations<-100 #number of iterations in the chain
	moveprobs<-c(1) # having length 1 disallows the new edge reversal move
	if(!(length(moveprobs)==1)){print('Vector of move probabilities has the wrong length!')}
  },
  "2"={ # with new edge reversal
	iterations<-100 #number of iterations in the chain
# Choose the probability of the different moves
# 1 is standard structure MCMC (including the possibility to stay still [officially needed for convergence])
# 2 is new edge reversal move
	moveprobs<-c(0.93,0.07)
	moveprobs<-moveprobs/sum(moveprobs) # normalisation
	if(!(length(moveprobs)==2)){print('Vector of move probabilities has the wrong length!')}
  },
  "3"={ # order MCMC
	iterations<-100 #number of iterations in the chain
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
  },
  "4"={ # partition MCMC
	iterations<-100 #number of iterations in the chain
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
  },
  "5"={ # partition MCMC with edge reversal
	iterations<-100 #number of iterations in the chain
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
  },
  {# if none is chosen, we have a problem
    print('Not implemented')
  })

stepsave<-1 #stepsave<-iterations/1000 #how often to save the result

  if(seedset>0){
    set.seed(seednumber) # choose one?
  }

# Run the MCMC codes

  starttime<-proc.time() # for timing the problem

  if(MCMCtype<3){

    startDAG<-matrix(numeric(n*n),nrow=n) # starting DAG is empty say
    revallowed<-1 # allow standard edge reversals

    example<-structureMCMC(n,startDAG,iterations,stepsave,maxparents,parenttable,scoretable,revallowed,moveprobs)

  } else if (MCMCtype==3){

    startorder<-c(1:n) # starting order

    example<-orderMCMC(n,startorder,iterations,stepsave,parenttable,scoretable,moveprobs)

  } else if (MCMCtype>3){

    startpermutation<-c(1:n) # pick a starting permutation
    startpartition<-c(n) # and a starting partition - c(n) gives the empty DAG

    example<-partitionMCMC(n,startpermutation,startpartition,iterations,stepsave,parenttable,scoretable,moveprobs)

  }

  endtime<-proc.time()
  endtime<-endtime-starttime
  print('Time to run the')
  print(iterations)
  print('iterations')
  print(endtime)

# Plot the results

  nparts<-length(example[[2]])
  DAGscores<-unlist(example[[2]])
  shiftedscores<-DAGscores
  maxscore<-max(shiftedscores)

  maxDAG<-example[[1]][which(DAGscores==max(DAGscores))][[1]]

# plotting parameters

  par(mar=c(2.25,4.75,1.25,0.5))
  par(mgp=c(3,1,0))
  par(cex.axis=1.25)
  par(cex.lab=1.5)

  plot(1:nparts,shiftedscores,type="l", ylab="DAG score", xlab="",#xlab="iteration step", 
main="", col="blue",ylim=c(maxscore-3.2,maxscore+0.2))
  abline(h=maxscore,col='springgreen4', lwd=3)
 # abline(h=0,col='red',lty=3,lwd=3)
  lines(1:nparts,shiftedscores,type="l", col="blue")

  if((saveoutput==1)&&(seedset==1)){
    switch(as.character(MCMCtype),
    "1"={ # standard structure MCMC
      typename<-"structurer"
      filetemp<-paste(dirname,typename,n,"nodesits",iterations/1000,"seed",seednumber,filetypes[filetypey],sep="")
	save(DAGscores,file=paste(dirname,"structure",n,"nodesscoresits",iterations/1000,"seed",seednumber,".RData",sep=""))
    },
    "2"={ # with new edge reversal
      typename<-"edgerev"
      filetemp<-paste(dirname,typename,n,"nodesone",100*moveprobs[1],"its",iterations/1000,"seed",seednumber,filetypes[filetypey],sep="")
	save(DAGscores,file=paste(dirname,"edgerev",n,"nodesscoresone",100*moveprobs[1],"its",iterations/1000,"seed",seednumber,".RData",sep=""))
    },
    "3"={ # order MCMC
  typename<-"order"
  filetemp<-paste(dirname,typename,n,"nodesone",100*moveprobs[1],"two",100*moveprobs[2],"its",iterations/1000,"seed",seednumber,filetypes[filetypey],sep="")  
	save(DAGscores,file=paste(dirname,"order",n,"nodesscoresone",100*moveprobs[1],"two",100*moveprobs[2],"its",iterations/1000,"seed",seednumber,".RData",sep=""))
    },
    "4"={ # partition MCMC
      typename<-"partition"
      filetemp<-paste(dirname,typename,n,"nodesone",100*moveprobs[1],"two",100*moveprobs[2],"three",100*moveprobs[3],"its",iterations/1000,"seed",seednumber,filetypes[filetypey],sep="")
	save(DAGscores,file=paste(dirname,"partition",n,"nodesscoresone",100*moveprobs[1],"two",100*moveprobs[2],"three",100*moveprobs[3],"its",iterations/1000,"seed",seednumber,".RData",sep=""))
    },
    "5"={ # partition MCMC with edge reversal
      typename<-"partedgerev"
      filetemp<-paste(dirname,typename,n,"nodesone",100*moveprobs[1],"two",100*moveprobs[2],"three",100*moveprobs[3],"six",100*moveprobs[6],"its",iterations/1000,"seed",seednumber,filetypes[filetypey],sep="")
	save(DAGscores,file=paste(dirname,"partition",n,"nodesedgerevscoresone",100*moveprobs[1],"two",100*moveprobs[2],"three",100*moveprobs[3],"six",100*moveprobs[6],"its",iterations/1000,"seed",seednumber,".RData",sep=""))
    },
    {# if none is chosen, we have a problem
	print('Not implemented!')
    })

  switch(as.character(filetypey),
   "1"={# png format
     png(filetemp, width=7.5, height=3.75, unit="in",res=100,bg="transparent")
   },
   "2"={# pdf format
     pdf(filetemp, width=7.5, height=3.75, onefile=F, pointsize=10,  paper="special")
   })
	
  par(mar=c(2.25,4.75,1.25,0.5))
  par(mgp=c(3,1,0))
  par(cex.axis=1.25)
  par(cex.lab=1.5)
	plot(1:nparts,shiftedscores,type="l", ylab="DAG score", xlab="",#xlab="iteration step", 
	     main="", col="blue",ylim=c(maxscore-3.2,maxscore+0.2))
	abline(h=maxscore,col='springgreen4', lwd=3)
	# abline(h=0,col='red',lty=3,lwd=3)
	lines(1:nparts,shiftedscores,type="l", col="blue")
	dev.off()
  
  }
 
createmovie<-0 # build a movie for presentations?

if(createmovie>0){

  imagemagickcommand<-"DAGmovie:\n\tconvert -delay 15 -dispose Background"
  for(i in 1:(iterations+1)){
    incidence<-example[[1]][i][[1]]
    filetemp<-"dag3nodes"
    for(j in 1:n){
      childnodes <- which(incidence[j,]==1)
      numchildren<-length(childnodes)
      filetemp<-paste(filetemp,paste(as.vector(rbind(rep(j,numchildren),childnodes)),collapse=""),sep="")
    }
    imagemagickcommand<-paste(imagemagickcommand," ",filetemp,".gif",sep="")
  }
  imagemagickcommand<-paste(imagemagickcommand," -loop 0 ",typename,"3nodemovie.gif",sep="")

  anothercommand<-"\n\ngraphmovie:\n\tconvert -delay 15 -dispose Background"
  for(i in 1:(iterations+1)){
    filetemp<-paste(dirname,typename,i,".png",sep="")  
    anothercommand<-paste(anothercommand," ",filetemp,sep="") 
    png(filetemp, width=7.5, height=3.75, unit="in",res=100,bg="transparent")
    par(mar=c(2.25,4.75,1.25,0.5))
    par(mgp=c(3,1,0))
    par(cex.axis=1.25)
    par(cex.lab=1.5)
    plot(1:nparts,shiftedscores,type="l", ylab="DAG score", xlab="",#xlab="iteration step", 
     main="", col="blue",ylim=c(maxscore-3.2,maxscore+0.2))
    abline(h=maxscore,col='springgreen4', lwd=3)
# abline(h=0,col='red',lty=3,lwd=3)
    lines(1:nparts,shiftedscores,type="l", col="blue")
    points(i,shiftedscores[i],col="red",pch=19,cex=2)
    dev.off()
  }
  anothercommand<-paste(anothercommand," -loop 0 ",typename,"graphmovie.gif",sep="")
  writeLines(paste(imagemagickcommand,anothercommand,sep=""),"makefile")
}