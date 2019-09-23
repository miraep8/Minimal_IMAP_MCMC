
# These are functions used by the structure code
# this code is taken from the Dortmund course programmed by Miriam Lohr

### calculation of the first ancestor matrix:
ancestor <- function(incidence){
incidence1 <- incidence
incidence2 <- incidence
k <- 1
while (k < nrow(incidence)){
incidence1 <- incidence1%*%incidence
incidence2 <- incidence2 + incidence1
k <-k+1
}
incidence2[which(incidence2[,]>0)] <- 1
return(t(incidence2))}

top_order <- function(incidence){
Order <- numeric(n)
fan_in <- numeric(n)
no_fan_in <- numeric(0)
m <- 1
for (p in 1:n){                                       # number of parent nodes at the beginning
fan_in[p] <- sum(incidence[,p])
}
no_fan_in <- which(fan_in==0)
while (length(which(Order==0))>0){                    # as long as there is a node without an order
fan_in[which(incidence[no_fan_in[1],]==1)] <- fan_in[which(incidence[no_fan_in[1],]==1)] - 1
no_fan_in <- c(no_fan_in, c(which(incidence[no_fan_in[1],]==1),which(fan_in==0))[duplicated(c(which(incidence[no_fan_in[1],]==1),which(fan_in==0)))])
Order[m] <- no_fan_in[1]
no_fan_in <- no_fan_in[-1]
m <- m+1
}
return(Order)
}

### assign the topological order of the descendants of the child
des_top_order <- function(incidence, ancest1,child){
top <- top_order(incidence)
position_child <- which(top==child)
top_all_after <- top[position_child:n]                # top. order without the "first" nodes
desc <- which(ancest1[,child]==1)                     # descendants of the child
inter_step <- c(child,desc,top_all_after)
des_top <- inter_step[which(duplicated(inter_step))]
return(des_top)
}

