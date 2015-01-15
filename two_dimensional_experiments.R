###########################################
####### Two-Dimensional Experiments #######
###########################################

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

N <- 1000
p <- 5
nu <- 0.05

temp_ladder <- exp(seq(0,log(nu),length=p+1))
prop_sd_ladder <- 1/temp_ladder
burnin <- 2000

runs <- 100

x_0 <- matrix(nrow=runs,ncol=2)
wpt_props <- matrix(nrow=runs,ncol=length(wts))
wpt_region <- vector('numeric',length=runs)

parallel_props <- matrix(nrow=runs,ncol=length(wts))
parallel_region <- vector('numeric',length=runs)

cl <- makeCluster(p)

for(i in 1:runs){
  x_0[i,] <- runif(2,-20,20)
  wpt_samps <- wpt(N=N,log_target=log_target,x_0=x_0[i,],p=p,nu=nu,prop_sd=prop_sd_ladder[p+1],burnin=burnin,parallel=TRUE,cl=cl,keep.open=TRUE)
  wpt_props[i,] <- colMeans(withinBounds(t(wpt_samps),grps))
  wpt_region[i] <- sum(wpt_props[i,])
  
  if(i%%(runs/10)==0) print(i)
}

stopCluster(cl)

for(i in 1:runs){
  
  parallel_samps <- parallel_tempering(N=N,log_target=log_target,x_0=x_0[i,],p=p,temp_ladder=temp_ladder,prop_sd_ladder=prop_sd_ladder,burnin=burnin)
  
  parallel_props[i,] <- colMeans(withinBounds(t(parallel_samps),grps))
  parallel_region[i] <- sum(parallel_props[i,])
  
  if(i%%(runs/10)==0) print(i)
  
}
