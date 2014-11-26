#######################################################
################## Truncated Normal ###################
#######################################################

rtruncnorm <- function(n,mean=0,sd=1,lower_bound=-Inf,upper_bound=Inf){
  
  rand <- runif(n,min=pnorm(lower_bound,mean=mean,sd=sd),max=pnorm(upper_bound,mean=mean,sd=sd))
  return(qnorm(rand,mean=mean,sd=sd))
  
}

dtruncnorm <- function(x,mean=0,sd=1,lower_bound=-Inf,upper_bound=Inf){
  
  return(dnorm(x,mean=mean,sd=sd)/(1-(pnorm(lower_bound,mean=mean,sd=sd)+pnorm(upper_bound,mean=mean,sd=sd,lower.tail=FALSE))))
  
}

ptruncnorm <- function(x,mean=0,sd=1,lower_bound=-Inf,upper_bound=Inf,lower.tail=TRUE){
  result <- vector('numeric',length(x))
  for(i in 1:length(x)){
    
    if(x[i]>=lower_bound&x[i]<=upper_bound){
      result[i] <- pnorm(x[i],mean=mean,sd=sd,lower.tail=lower.tail)+pnorm(upper_bound*lower.tail+lower_bound*!lower.tail,mean=mean,sd=sd,lower.tail=!lower.tail)
    } else if(x[i]>upper_bound){
      result[i] <- lower.tail
    } else {
      result[i] <- !lower.tail
    }
    
  }
  return(result)
}