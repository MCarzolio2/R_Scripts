###########################################
###### KS Convergence Experiments #########
###########################################

source('./weighted_particle_tempering.R')
source('./parallel_tempering.R')

log_target <- function(x){
  mus <- c(0,10)
  sds <- c(1,2)
  wts <- c(0.2,0.8)
  log_terms <- -1/2*(log(2*pi)+2*log(sds)+(x-mus)^2/sds^2)
  return(log(wts%*%exp(log_terms)))
  
}

target_cdf <- function(x){
  mus <- c(0,10)
  sds <- c(1,2)
  wts <- c(0.2,0.8)
  probs <- sapply(x,pnorm,mean=mus,sd=sds)
  return(wts%*%probs)
}

set.seed(123)
p <- 8
N <- 1000
runs <- 100

# WPT parameters
nu <- 0.01
prop_sd <- 1
delta <- nu

# Parallel Tempering parameters
temp_ladder <- exp(seq(0,log(nu),length=p+1))
prop_sd_ladder <- seq(0.1,5,length=p+1)

# Initializations for algorithms
x_0 <- runif(runs,-20,20)

#setwd('C:/Users/mcarzolio/Documents/Projects/Ultra Tempering Algorithm/Experiments')
wpt_filename <- 'wpt_samples.csv'
pt_filename <- 'pt_samples.csv'

# Initialize
cl <- makeCluster(p)
varlist <- c('N','log_target','x_0','prop_sd','rwmh','nu')
clusterExport(cl,varlist)

write.table(wpt(N=N,log_target=log_target,x_0=x_0[1],p=p,nu=nu,prop_sd=prop_sd,delta=delta,burnin=0,parallel=TRUE,cl=cl,keep.open=TRUE),wpt_filename,row.names=FALSE,col.names=FALSE,sep=',')
write.table(parallel_tempering(N=N,log_target=log_target,x_0=x_0[1],p=p,temp_ladder=temp_ladder,prop_sd_ladder=prop_sd_ladder,burnin=0),pt_filename,row.names=FALSE,col.names=FALSE,sep=',')

# Loop through all runs
for(i in 2:runs){
  write.table(wpt(N=N,log_target=log_target,x_0=x_0[i],p=p,nu=nu,prop_sd=prop_sd,delta=delta,burnin=0,parallel=TRUE,cl=cl,keep.open=TRUE),wpt_filename,row.names=FALSE,col.names=FALSE,sep=',',append=TRUE)
  write.table(parallel_tempering(N=N,log_target=log_target,x_0=x_0[i],p=p,temp_ladder=temp_ladder,prop_sd_ladder=prop_sd_ladder,burnin=0),pt_filename,row.names=FALSE,col.names=FALSE,sep=',',append=TRUE)
  if(i%%(runs/10)==0) print(i)
}

stopCluster(cl)

# Import all samples
wpt_samples <- read.csv(wpt_filename,header=FALSE)
pt_samples <- read.csv(pt_filename,header=FALSE)

# Compute KS statistics
wpt_ks_stats <- vector('numeric',N)
pt_ks_stats <- vector('numeric',N)

for(i in 1:N){
  
  wpt_ks_stats[i] <- ks.test(wpt_samples[,i],target_cdf)$statistic
  pt_ks_stats[i] <- ks.test(pt_samples[,i],target_cdf)$statistic
  
}

# Produce plots
plot(smooth.spline(wpt_ks_stats),type='l',lwd=2,ylim=c(0,0.4),xlab='Iteration',ylab='Smoothed Spline Over KS-Statistic',main=substitute(paste(delta, ' = ' ,i), list(i = delta)))
lines(smooth.spline(pt_ks_stats),col='black',lwd=2,lty=2)
