###########################################
########### Parallel Tempering ############
###########################################

source('./random_walk_metropolis_hastings.R')

parallel_tempering <- function(N,log_target,x_0,p,temp_ladder,prop_sd_ladder,burnin=0,high_dim_exp=FALSE){
  
  d <- length(x_0)
  
  samples <- matrix(ncol=N,nrow=d)
  
  current_sample <- if(!high_dim_exp) matrix(x_0,ncol=p+1,nrow=d) else matrix(runif((p+1)*d,0,10),ncol=p+1,nrow=d)
  
  temp_ladder <- if(1%in%temp_ladder) sort(temp_ladder,decreasing=TRUE) else c(1,sort(temp_ladder,decreasing=TRUE))
  #swap_rate <- vector('logical',p)
  #times_sampled <- vector('logical',p)
  for(i in 1:(N-1+burnin)){
    if(i>burnin){
      samples[,i-burnin] <- current_sample[,1]
    }
    
    # Swap randomly selected pair of particles:
    random_ind <- sample(1:p,1)
    #times_sampled[random_ind] <- times_sampled[random_ind]+1
    log_mh_ratio <- (temp_ladder[random_ind+1]-temp_ladder[random_ind])*(log_target(current_sample[,random_ind])-log_target(current_sample[,random_ind+1]))
    if(log(runif(1))<log_mh_ratio){
      current_sample[,random_ind:(random_ind+1)] <- current_sample[,(random_ind+1):random_ind]
      #swap_rate[random_ind] <- swap_rate[random_ind]+1
    }
    
    # Evolve other particles:
    for(j in (1:(p+1))[-c(random_ind:(random_ind+1))]){
      
      current_sample[,j] <- rwmh(N = 1,log_target = function(x) temp_ladder[j]*log_target(x),x_0 = current_sample[,j],prop_sd = prop_sd_ladder[j],burnin = 0)
      
    }
    
  }
  samples[,N] <- current_sample[,1]
  return(samples)
  
}