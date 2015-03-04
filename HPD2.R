######### Highest Posterior Density
if(!'mclust'%in%installed.packages()[,1]) install.packages('mclust')
library(mclust)
getIndex <- function(base,num,len){
	
	result <- if(floor(log(num,base))>-Inf) vector('numeric',len) else return(rep(0,len))
	temp <- num
	while(temp>0){
		result[floor(log(temp,base))+1] <- floor(temp/base^floor(log(temp,base)))
		temp <- temp - floor(temp/base^floor(log(temp,base)))*base^(floor(log(temp,base)))
		}
	return(result)
	
	}


HPD <- function(bottom,top,width,func){
	
	# Takes a bottom-most corner ('bottom') and top-most corner ('top') of the p-dimensional support of 'func'
	# chops the resulting hypercube into bins of width 'width'
	# returns a dataframe where each row corresponds to a bin
	# the first column of the dataframe is an index such that, when written out in base-p, each digit refers to a coordinate in p-space, counting the number of bins in each direction one must travel to find that particular row's bin
	# the next p columns specify the coordinate of the bottom-most corner of the bin
	# the second to last column is the approximate measure of that bin
	# the last column is the cumulative approximate measure of all bins at and above the given row
	p <- length(bottom)
	numBins <- sapply((top-bottom)/width,ceiling)
	result <- matrix(nrow=prod(numBins),ncol=p+3)
	result[,1] <- 1:nrow(result) - 1
	base <- max(numBins)
	
	for(i in 1:nrow(result)){
		result[i,2:(p+1)] <- bottom + getIndex(base,result[i,1],p)*width
		result[i,p+2] <- func(result[i,2:(p+1)])
		}
	result <- result[order(result[,p+2],decreasing=T),]	
	result[,p+3] <- cumsum(result[,p+2])
	return(data.frame(result))
	}
	
consolidateGroup <- function(grp,width){
	
	newGrp <- rbind(t(apply(grp,1,FUN=function(x) x+c(width/2,width/2))),t(apply(grp,1,FUN=function(x) x-c(width/2,width/2))),t(apply(grp,1,FUN=function(x) x+c(-width/2,width/2))),t(apply(grp,1,FUN=function(x) x+c(width/2,-width/2))))
	
	inds1 <- which(duplicated(newGrp))
	inds2 <- which(duplicated(newGrp,fromLast=T))
	inds <- unique(c(inds1,inds2))
	if(length(inds)>0) return(newGrp[-inds,]) else return(grp)
	
	}	
	
groups <- function(bottoms,width){
	# determines which bins, defined by 'bottoms' and 'width' are grouped together
	# now we are only working in 2D
	
	#clust <- Mclust(bottoms)
	clust <- kmeans(bottoms,centers=t(matrix(c(-5,-8,5,5,10,12,-15,5,5,-15), nrow=2,byrow=F)))
	print(clust)
	grps <- list(NULL)
	for(i in 1:5){ #clust$G){
		#grps <- c(grps,list(bottoms[clust$classification==i,]))
		grps <- c(grps,list(bottoms[clust$cluster==i,]))
		grps[[i+1]] <- consolidateGroup(grps[[i+1]],width)
		}	
	
	return(grps[-1])
	}	
	
hpd  <- HPD(bottom=c(-20,-20),top=c(20,20),width=0.1,target)
	
	
plotHPD <- function(hpd,plane=1:2,width,probs,ret=T,grps=NULL){
	if(is.null(grps)){
			probs <- c(0,sort(probs))
		for(i in 2:length(probs)){
			prob1 <- probs[i]*hpd[nrow(hpd),ncol(hpd)]
			prob2 <- probs[i-1]*hpd[nrow(hpd),ncol(hpd)]
			set <- hpd[hpd[,ncol(hpd)]<=prob1&hpd[,ncol(hpd)]>prob2,1+plane]
			grps <- groups(set,width)
	
			}
			
		}
	for(i in 1:length(grps)){
	  
	  temp <- data.frame(grps[[i]])
	  x <- mean(temp[,1])
	  y <- mean(temp[,2])
	  
	  angs <- atan((temp[,2]-y)/(temp[,1]-x))
	  angs1 <- angs[temp[,1]>=x&temp[,2]>=y]
	  angs2 <- angs[temp[,1]<x&temp[,2]>=y]
	  angs3 <- angs[temp[,1]<x&temp[,2]<y]
	  angs4 <- angs[temp[,1]>=x&temp[,2]<y]
	  
	  temp1 <- temp[temp[,1]>=x&temp[,2]>=y,][order(angs1),]
	  temp2 <- temp[temp[,1]<x&temp[,2]>=y,][order(angs2),]
	  temp3 <- temp[temp[,1]<x&temp[,2]<y,][order(angs3),]
	  temp4 <- temp[temp[,1]>=x&temp[,2]<y,][order(angs4),]
	  
	  newtemp <- rbind(temp1,temp2,temp3,temp4)
	  newtemp <- rbind(newtemp,newtemp[1,])
	  
	  lines(newtemp,lwd=2,lty=2)
	  grps[[i]] <- newtemp
	}
	if(ret) return(grps)
}
plot(0,type='n')
plotHPD(hpd,width=0.1,probs=0.9)->grps

withinBounds <- function(pts,bounds){
	# bounds needs to be a list of n x 2 matrices
	# each matrix in bounds defines the vertices of the boundaries of non-overlapping polygons
	precision <- 6
	ints <- list(NULL)
	for(i in 1:length(bounds)){
		ints <- c(ints,list(matrix(nrow=nrow(pts),ncol=nrow(bounds[[i]])-1)))
		}
	ints <- ints[-1]	
	
	for(i in 1:length(bounds)){
		
		for(j in 2:nrow(bounds[[i]])){
				coeffs <- lm(bounds[[i]][c(j-1,j),2]~bounds[[i]][c(j-1,j),1])$coefficients
				intY <- if(!is.na(coeffs[2])) round(coeffs[1]+coeffs[2]*pts[,1],precision) else Inf
				ints[[i]][intY<Inf,j-1] <- (round(intY-max(bounds[[i]][c(j-1,j),2]),precision)<=0)&(round(intY-min(bounds[[i]][c(j-1,j),2]),precision)>=0)&(round(pts[,1]-max(bounds[[i]][c(j-1,j),1]),precision)<=0)&(round(pts[,1]-min(bounds[[i]][c(j-1,j),1]),precision)>=0)&(round(pts[,1]-max(bounds[[i]][c(j-1,j),1]),precision)<=0)&(round(pts[,2]-intY,precision)<=0)
				
				ints[[i]][intY==Inf,j-1] <- 0
			
			}
		}	
	
	result <- sapply(1:length(ints),FUN=function(x) rowSums(ints[[x]])%%2==1)
	
	return(result)
	}	

trueProbs <- withinBounds(hpd[,2:3],grps)
w <- sapply(1:ncol(trueProbs),FUN=function(x) sum(hpd[trueProbs[,x],4])/hpd[nrow(hpd),ncol(hpd)])

if(interactive()){
  #system.time(landedParallel <- withinBounds(Zp,grps))
  #mean(rowSums(landedParallel))
  #cbind(sort(colSums(landedParallel)/nrow(landedParallel)/0.9),sort(w))
  #system.time(landedUltra <- withinBounds(Zu,grps))
  #mean(rowSums(landedUltra))
  #cbind(sort(colSums(landedUltra)/nrow(landedUltra)),sort(w))
  #plot(sapply(1:nrow(Zu),FUN=function(x)mean(landedUltra[1:x,5])))
}


