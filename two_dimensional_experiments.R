###########################################
####### Two-Dimensional Experiments #######
###########################################

#setwd("C:/Users/mcarzolio/Documents/GitHub/R_Scripts/")
source('./weighted_particle_tempering.R')
source('./parallel_tempering.R')

mus <- matrix(c(c(-5,-8),c(5,5),c(10,12),c(-15,5),c(5,-15)),nrow=2)
sds <- 1
wts <- c(1/2,1/6,1/12,1/6,1/12)

log_target <- function(x){
  mus <- matrix(c(c(-5,-8),c(5,5),c(10,12),c(-15,5),c(5,-15)),nrow=2)
  sds <- 1
  wts <- c(1/2,1/6,1/12,1/6,1/12)
  log_terms <- colSums(-1/2*(log(2*pi)+2*log(sds)+(x-mus)^2/sds^2))
  return(log(wts%*%exp(log_terms)))
}

target <- function(x){
  return(exp(log_target(x)))
}

source('./HPD2.R')
set.seed(123)
N <- 1000
burnin <- 2000
runs <- 100

num_particles <- c(2,5,10,20)
temps <- c(0.9,0.5,0.1,0.05,0.01,0.005)

wpt_mse <- matrix(0,ncol=length(temps),length(num_particles))
pt_mse <- matrix(0,ncol=length(temps),length(num_particles))

for(nu in temps){
  for(p in num_particles){
    
    temp_ladder <- exp(seq(0,log(nu),length=p+1))
    prop_sd_ladder <- temp_ladder^(-0.5)
    
    x_0 <- matrix(nrow=runs,ncol=2)
    wpt_props <- matrix(nrow=runs,ncol=length(wts))
    wpt_region <- vector('numeric',length=runs)
    
    parallel_props <- matrix(nrow=runs,ncol=length(wts))
    parallel_region <- vector('numeric',length=runs)
    
    cl <- makeCluster(p)
    
    for(i in 1:runs){
      x_0[i,] <- runif(2,-20,20)
      wpt_samps <- wpt(N=N,log_target=log_target,x_0=x_0[i,],p=p,nu=nu,prop_sd=prop_sd_ladder[p+1],delta=nu,burnin=burnin,parallel=TRUE,cl=cl,keep.open=TRUE)
      wpt_props[i,] <- colMeans(withinBounds(t(wpt_samps),grps))
      wpt_region[i] <- sum(wpt_props[i,])
      #wpt_mse[which(num_particles==p),which(temps==nu)] <- wpt_mse[which(num_particles==p),which(temps==nu)]+sum((wpt_props[i,]-w)^2)
      if(i%%(runs/10)==0) print(paste('WPT iteration ',i,' for nu = ',nu,' and p = ',p,sep=''))
    }
  
    stopCluster(cl)
    
    #wpt_mse[which(num_particles==p),which(temps==nu)] <- wpt_mse[which(num_particles==p),which(temps==nu)]/runs
    
    for(i in 1:runs){
      
      parallel_samps <- parallel_tempering(N=N,log_target=log_target,x_0=x_0[i,],p=p,temp_ladder=temp_ladder,prop_sd_ladder=prop_sd_ladder,burnin=burnin)
      
      parallel_props[i,] <- colMeans(withinBounds(t(parallel_samps),grps))
      parallel_region[i] <- sum(parallel_props[i,])
      #pt_mse[which(num_particles==p),which(temps==nu)] <- pt_mse[which(num_particles==p),which(temps==nu)]+sum((parallel_props[i,]-w)^2)
      if(i%%(runs/10)==0) print(paste('PT iteration ',i,' for nu = ',nu,' and p = ',p,sep=''))
      
    }
    #pt_mse[which(num_particles==p),which(temps==nu)] <- pt_mse[which(num_particles==p),which(temps==nu)]/runs
    
    if(p==num_particles[2]&nu==temps[length(temps)-1]){
      keep_parallel_props <- parallel_props
      keep_parallel_region <- parallel_region
    }
    if(p==num_particles[2]&nu==temps[length(temps)-2]){
      keep_wpt_props <- wpt_props
      keep_wpt_region <- wpt_region
    }
  }
}
