###########################################
####### High-Dimensional Experiment #######
###########################################
setwd('~/GitHub/R_Scripts')
source('./weighted_particle_tempering.R')
source('./parallel_tempering.R')

log_target <- function(x){
  
  d <- length(x)
  exp_term1 <- -1/2*t(x)%*%x
  exp_term2 <- -1/2*t(x-rep(10,d))%*%(x-rep(10,d))
  max_term <- max(exp_term1,exp_term2)
  
  result <- max_term + log(exp(exp_term1-max_term)+exp(exp_term2-max_term))
  return(result)
  
}

#############
#  dim = 5: #
#############
set.seed(123)
N <- 1000
d <- 5
x_0 <- runif(d,0,10)

# Tune RWMH first:
nu <- 0.1
prop_sd <- 3
rwmh_samps <- rwmh(N,function(x) nu*log_target(x),x_0,prop_sd)
matplot(t(rwmh_samps),type='l')

# Find the lowest p that yields good results:
p <- 10
wpt_samps <- wpt(N,log_target,x_0,p,nu,prop_sd,delta=nu,high_dim_exp=TRUE)[[1]]
matplot(t(wpt_samps),type='l')

# Now try parallel tempering:
pt_samps <- parallel_tempering(N,log_target,x_0,p,temp_ladder=exp(seq(0,log(nu),length=p+1)),prop_sd_ladder=exp(-seq(0,log(nu),length=p+1)))
matplot(t(pt_samps),type='l')

#############
# dim = 50: #
#############

set.seed(123)
N <- 1000
d <- 50
x_0 <- runif(d,0,10)

# Tune RWMH first:
nu <- 0.5
prop_sd <- 0.5
rwmh_samps <- rwmh(N,function(x) nu*log_target(x),x_0,prop_sd)
matplot(t(rwmh_samps),type='l',xlab='Iteration',main='50-Dimensional RWM Samples',ylab='',lwd=2,sub=expression((nu==0.5)))

# Find the lowest p that yields good results:
p <- 50
wpt_samps <- wpt(N,log_target,x_0,p,nu,prop_sd,delta=nu,high_dim_exp=TRUE)[[1]]
matplot(t(wpt_samps),type='l',xlab='Iteration',main='50-Dimensional WPT Samples',ylab='',lwd=2,sub=expression((nu==0.5)))

# Tune RWMH first:
nu <- 0.9
prop_sd <- 0.5
rwmh_samps <- rwmh(N,function(x) nu*log_target(x),x_0,prop_sd)
matplot(t(rwmh_samps),type='l',xlab='Iteration',main='50-Dimensional RWM Samples',ylab='',lwd=2,sub=expression((nu==0.03)))

# Find the lowest p that yields good results:
p <- 50
wpt_samps <- wpt(N,log_target,x_0,p,nu,prop_sd,delta=0,high_dim_exp=TRUE)[[1]]
matplot(t(wpt_samps),type='l',xlab='Iteration',main='50-Dimensional WPT Samples',ylab='',lwd=2,sub=expression((nu==0.03)))

# Now try parallel tempering:
p <- 20
nu <- 0.01
pt_samps <- parallel_tempering(N,log_target,x_0,p,temp_ladder=exp(seq(0,log(nu),length=p+1)),prop_sd_ladder=exp(-seq(0,log(nu),length=p+1))/5,high_dim_exp=TRUE)
matplot(t(pt_samps),type='l',xlab='Iteration',main='50-Dimensional PT Samples',ylab='',lwd=2)


#############
# dim = 100: #
#############

set.seed(123)
N <- 10000
d <- 100
x_0 <- runif(d,0,10)

# Tune RWMH first:
nu <- 0.9
prop_sd <- 0.5
rwmh_samps <- rwmh(N,function(x) nu*log_target(x),x_0,prop_sd)
matplot(t(rwmh_samps),type='l',xlab='Iteration',main='50-Dimensional RWM Samples',ylab='',lwd=2,sub=expression((nu==0.03)))

# Find the lowest p that yields good results:
p <- 10
wpt_samps <- wpt(N,log_target,x_0,p,nu,prop_sd,delta=0,high_dim_exp=TRUE)
matplot(t(wpt_samps[[1]]),col=rgb(190/255,190/255,190/255,0.1),type='l',ylab='',xlab='Iteration',main='100-Dimensional WPT Samples')
mtext(substitute(paste("(",nu==n,", ",p==pee,")"),list(n=nu,pee=p)))
points(wpt_samps[[1]][1,],col=rgb(1,0,0,0.2),pch=16)
