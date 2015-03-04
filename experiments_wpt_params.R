###########################################
###### Experimenting with parameters ######
###########################################


source('./weighted_particle_tempering.R')
log_target <- function(x){
  
  d <- length(x)
  exp_term1 <- -1/2*t(x)%*%x
  exp_term2 <- -1/2*t(x-rep(10,d))%*%(x-rep(10,d))
  max_term <- max(exp_term1,exp_term2)
  
  result <- max_term + log(exp(exp_term1-max_term)+exp(exp_term2-max_term))
  return(result)
  
}

###############
# delta = nu: #
###############

set.seed(123)
N <- 1000
d <- 10
x_0 <- runif(d,0,10)
p <- c(20,50,100,200,500)

lmhr <- matrix(nrow=N,ncol=length(p))
delta <- nu
for(j in 1:length(p)){
  wpt_samps <- wpt(N,log_target,x_0,p[j],nu,prop_sd,delta,high_dim_exp=TRUE,parallel=FALSE)
  lmhr[,j] <- wpt_samps[[2]]
}
boxplot(lmhr,axes=FALSE,xlab='Number of particles',ylab='Log MH Ratio')
axis(1,label=p,at=1:length(p))
axis(2)

###############
# delta = 1:  #
###############
# Corresponds to picking particles uniformly
lmhr <- matrix(nrow=N,ncol=length(p))
delta <- 1
for(j in 1:length(p)){
  wpt_samps <- wpt(N,log_target,x_0,p[j],nu,prop_sd,delta,high_dim_exp=TRUE,parallel=FALSE)
  lmhr[,j] <- wpt_samps[[2]]
}
boxplot(lmhr,axes=FALSE,xlab='Number of particles',ylab='Log MH Ratio')
axis(1,label=p,at=1:length(p))
axis(2)

###############
# delta = 0:  #
###############
# Corresponds to picking particles with weights prop to target dist
lmhr <- matrix(nrow=N,ncol=length(p))
delta <- nu
for(j in 1:length(p)){
  wpt_samps <- wpt(N,log_target,x_0,p[j],nu,prop_sd,delta,high_dim_exp=TRUE,parallel=FALSE)
  lmhr[,j] <- wpt_samps[[2]]
}
boxplot(lmhr,axes=FALSE,xlab='Number of particles',ylab='Log MH Ratio')
axis(1,label=p,at=1:length(p))
axis(2)

###############
# delta = .5: #
###############
# Corresponds to picking particles with weights prop to target dist
lmhr <- matrix(nrow=N,ncol=length(p))
delta <- 0.5
for(j in 1:length(p)){
  wpt_samps <- wpt(N,log_target,x_0,p[j],nu,prop_sd,delta,high_dim_exp=TRUE,parallel=FALSE)
  lmhr[,j] <- wpt_samps[[2]]
}
boxplot(lmhr,axes=FALSE,xlab='Number of particles',ylab='Log MH Ratio')
axis(1,label=p,at=1:length(p))
axis(2)

###############
# delta = .9: #
###############
# Corresponds to picking particles with weights prop to target dist
lmhr <- matrix(nrow=N,ncol=length(p))
delta <- 0.9
for(j in 1:length(p)){
  wpt_samps <- wpt(N,log_target,x_0,p[j],nu,prop_sd,delta,high_dim_exp=TRUE,parallel=FALSE)
  lmhr[,j] <- wpt_samps[[2]]
}
boxplot(lmhr,axes=FALSE,xlab='Number of particles',ylab='Log MH Ratio')
axis(1,label=p,at=1:length(p))
axis(2)

