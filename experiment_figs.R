###########################################
####### Produce Plots for WPT Paper #######
###########################################

# Two dimensional experiments:
#source('./two_dimensional_experiments.R')
fig_directory <- 'C:/Users/mcarzolio/Documents/Projects/Ultra Tempering Algorithm/Plots from Experiments/'
width <- 600
height <- 600

png(paste(fig_directory,'90_percent_coverage.png',sep=''),width,height)
boxplot(keep_parallel_region,keep_wpt_region,notch=TRUE,lwd=2,names=c('Parallel Tempering','Weighted Particle Tempering'))
abline(h=sum(w),lwd=2,lty=2)
dev.off()

png(paste(fig_directory,'wpt_coverage_props.png',sep=''),width,height)
boxplot(keep_wpt_props[,order(w)],notch=TRUE,lwd=2,names=sort(round(w,2)),ylim=c(0,1))
for(i in 1:length(w)) lines(c(0.5+i-1,1+i-0.5),c(sort(w)[i],sort(w)[i]),lwd=2,lty=2)
dev.off()

png(paste(fig_directory,'parallel_coverage_props.png',sep=''),width,height)
boxplot(keep_parallel_props[,order(w)],notch=TRUE,lwd=2,names=sort(round(w,2)),ylim=c(0,1))
for(i in 1:length(w)) lines(c(0.5+i-1,1+i-0.5),c(sort(w)[i],sort(w)[i]),lwd=2,lty=2)
dev.off()

if(!require(xtable)){
  install.packages('xtable')
  library(xtable)
}
xtable(data.frame(sqrt(pt_mse)),digits=3)
xtable(data.frame(sqrt(wpt_mse)),digits=3)

png(paste(fig_directory,'parallel_samps.png',sep=''),width,height)
plot(t(parallel_samps),pch='.',xlab='',ylab='',cex=3,ylim=c(-20,20),xlim=c(-20,20))
plotHPD(hpd,width=0.1,probs=0.9,grps=grps,ret=FALSE)
dev.off)

png(paste(fig_directory,'wpt_samps.png',sep=''),width,height)
plot(t(wpt_samps),pch='.',xlab='',ylab='',cex=3,ylim=c(-20,20),xlim=c(-20,20))
plotHPD(hpd,width=0.1,probs=0.9,grps=grps,ret=FALSE)
dev.off()

source('./ks_experiments.R')

png(paste(fig_directory,'ks_convergence.png',sep=''),width,height)
plot(smooth.spline(wpt_ks_stats),type='l',lwd=2,xlab='Iteration',ylab='Smoothed-Spline Over KS Statistic',ylim=c(0,0.4))
lines(smooth.spline(pt_ks_stats),lwd=2,lty=2)
dev.off()

x <- seq(-3,15,0.01)
plot(x,sapply(x,FUN=function(i) exp(log_target(i))),type='l',lwd=2,ylab=expression(pi(x)),xlab='x')
