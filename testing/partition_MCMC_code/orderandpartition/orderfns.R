# This function gives the list of scores compatible with each order

orderscore<-function(n,scorenodes,parenttable,scoretable,permy){

  orderscores<-rep(0,n)
  allscores<-vector("list",n)
  allowedscorerows<-vector("list",n)

  tablesize<-dim(parenttable[[1]]) # just to remove some arguments

  for (i in scorenodes){
    position<-which(permy==i)
    if(position==n){ # no parents are allowed
	orderscores[i]<-scoretable[[i]][1,1]
	allscores[[i]]<-orderscores[i] # there is only one score
	allowedscorerows[[i]]<-c(1) # there is only one score
    } else {
      if(position>1){
	bannednodes<-permy[c(1:(position-1))]
        allowedrows<-c(2:tablesize[1])
        for (j in 1:tablesize[2]){ # working columnwise allows R to speed up 
          bannedrows<-which(parenttable[[i]][allowedrows,j]%in%bannednodes)
	  if(length(bannedrows)>0){
	    allowedrows<-allowedrows[-bannedrows]
	  }
        }
	allowedrows<-c(1,allowedrows) # this row is always allowed
      	allscores[[i]]<-scoretable[[i]][allowedrows,1]
	allowedscorerows[[i]]<-allowedrows
      } else{ # all parents are allowed
	allscores[[i]]<-scoretable[[i]][,1]
	allowedscorerows[[i]]<-c(1:tablesize[1])
      }
      maxallowed<-max(allscores[[i]])
      orderscores[i]<-maxallowed+log(sum(exp(allscores[[i]]-maxallowed)))
    }
  }

  scores<-list()
  scores$allscores<-allscores
  scores$allowedrows<-allowedscorerows
  scores$totscores<-orderscores
return(scores)
}

