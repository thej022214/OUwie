#Does contour plot for likelihood surface for pair of parameters
library(akima)
library(lattice)
library(grDevices)
#written by Brian C. OMeara

OUwie.contour<-function(phy,data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"), simmap.tree=FALSE, scaleHeight=FALSE, root.station=TRUE, lb=0.000001, ub=1000, focal.param=NULL, clade=NULL, mserr="none", nrep=1000, sd.mult=3, levels=c(0.5,1,1.5,2),likelihood.boundary=Inf,lwd=2, ...){
#focal.param is something like c("alpha_2","sigma.sq_1"). They are then split on "_"
	if(length(focal.param)!=2) {
		stop("need a focal.param vector of length two")
	}
	if(sum(grepl("theta",focal.param))>0) {
		stop("contour mapping currently only works for alpha and sigma.sq parameters") 
	}
	globalMLE<-OUwie(phy=phy,data=data,model=model,simmap.tree=simmap.tree,scaleHeight=scaleHeight,root.station=root.station,lb=lb, ub=ub, clade=clade, mserr=mserr, diagn=TRUE)
	focal.param.df<-data.frame(strsplit(focal.param,"_"),stringsAsFactors=FALSE)
	names(focal.param.df)<-c(1,2)
	focal.param.df<-rbind(focal.param.df,rep(NA,2))
	focal.param.df<-rbind(focal.param.df,rep(NA,2))
	row.names(focal.param.df)<-c("parameter","element","MLE","SE")
	for (i in 1:2) {
		focal.param.df[3,i]<-as.numeric(globalMLE$solution[which(row.names(globalMLE$solution)==focal.param.df[1,i]),focal.param.df[2,i]])
		focal.param.df[4,i]<-as.numeric(globalMLE$solution.se[which(row.names(globalMLE$solution.se)==focal.param.df[1,i]),focal.param.df[2,i]])
	}
	rnorm.bounded<-function(mean=0,sd=1,bound=0) {
		x<-(-Inf)
		while (x<=bound) {
			x<-rnorm(1,mean,sd) 
		}
		return(x)
	}
	
	#now get our random points, sampling most densely near the MLE (likely where our contour lines will be drawn)
	param1.points<-replicate(n=round(nrep/4),expr=rnorm.bounded(mean=as.numeric(focal.param.df[3,1]),sd=as.numeric(focal.param.df[4,1])))
	param2.points<-replicate(n=round(nrep/4),expr=rnorm.bounded(mean=as.numeric(focal.param.df[3,2]),sd=as.numeric(focal.param.df[4,2])))
	param1.points<-c(param1.points,replicate(n=round(nrep/4),expr=rnorm.bounded(mean=as.numeric(focal.param.df[3,1]),sd=sd.mult*as.numeric(focal.param.df[4,1]))))
	param2.points<-c(param2.points,replicate(n=round(nrep/4),expr=rnorm.bounded(mean=as.numeric(focal.param.df[3,2]),sd=sd.mult*as.numeric(focal.param.df[4,2]))))
	
	#now lets figure out what the overall boundaries will be
	xlim=range(param1.points,as.numeric(focal.param.df[3,1])-1.96*as.numeric(focal.param.df[4,1]),as.numeric(focal.param.df[3,1])+1.96*as.numeric(focal.param.df[4,1]))
	ylim=range(param2.points,as.numeric(focal.param.df[3,2])-1.96*as.numeric(focal.param.df[4,2]),as.numeric(focal.param.df[3,2])+1.96*as.numeric(focal.param.df[4,2]))
	#I added this in there so that the hinterland does not include values that are not biologically possible -- JMB:
	if(xlim[1]<0){
		xlim[1]=0
	}
	if(ylim[1]<0){
		ylim[1]=0
	}
	if(strsplit(focal.param,"_")[[1]][1] == strsplit(focal.param,"_")[[2]][1]) {
		xlim=range(c(param1.points,param2.points,xlim,ylim))
		ylim=range(c(param1.points,param2.points,ylim,xlim))
	}
	
	#and sample the rest of the points from this space
	param1.points<-c(param1.points,runif(nrep-round(nrep/2),min=min(xlim),max=max(xlim)))
	param2.points<-c(param2.points,runif(nrep-round(nrep/2),min=min(ylim),max=max(ylim)))
	
	#now randomize order, in case we later want to do incremental saving or drawing
	param1.points<-sample(param1.points,size=length(param1.points),replace=FALSE)
	param2.points<-sample(param2.points,size=length(param2.points),replace=FALSE)
	
	dev<-function(p,phy,data,model,simmap.tree,scaleHeight,root.station,focal.param.vector,clade,mserr,globalMLE){ 
	#globalMLE is just for figuring out the structure of alpha and sigma.sq
		nRegimes<-dim(globalMLE$index.mat)[2]
		alpha<-rep(NA,nRegimes)
		sigma.sq<-rep(NA,nRegimes)
		for (freeParam in sequence(max(globalMLE$index.mat))) {
			entries<-which(globalMLE$index.mat==freeParam,arr.ind=TRUE)
			paramValue<-NA
			firstEntry<-entries[1,]
			paramNameRoot<-row.names(entries)[1]
			firstEntryName<-paste(paramNameRoot,firstEntry[2],sep="_")
			matchingSetParams<-which(firstEntryName==names(focal.param.vector))
			if (length(matchingSetParams)==1) {
				paramValue<-as.numeric(c(focal.param.vector[matchingSetParams]))
			}
			else {
				paramValue<-p[1]
				p<-p[-1] #pop off value we just used
			}
			for(i in sequence(dim(entries)[1])) {
				if(row.names(entries)[i]=="alpha") {
					alpha[entries[i,2]]<-paramValue
				}
				else {
					sigma.sq[entries[i,2]]<-paramValue
				}
			}
		}
		loglik=Inf
		try(loglik<-OUwie.fixed(phy=phy,data=data, model=model,simmap.tree=simmap.tree,scaleHeight=scaleHeight,root.station=root.station,alpha=alpha, sigma.sq=sigma.sq, theta=NULL, clade=clade, mserr=mserr, diagn=TRUE)$loglik)
		if(loglik==Inf){
			loglik=-10000000
		}
		print(paste("loglik is ",loglik))
		return(-loglik)
	}
	
	optimizeSemifixed<-function(X,phy,data,model,simmap.tree,scaleHeight,root.station,lb,ub,clade,mserr,globalMLE) {
		focal.param.vector<-X
		np<-max(globalMLE$index.mat)-length(focal.param.vector)
		lower = rep(lb, np)
		upper = rep(ub, np)
		ip<-1
		opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5, "xtol_rel"=.Machine$double.eps^0.5)
		out = nloptr(x0=rep(ip, length.out = np), eval_f=dev, opts=opts, phy=phy,data=data, model=model,simmap.tree=simmap.tree,scaleHeight=scaleHeight,root.station=root.station, lb=lower, ub=upper, clade=clade, mserr=mserr, globalMLE=globalMLE,focal.param.vector=focal.param.vector)
		return(-1*out$objective)
	}
	
	params.points<-data.frame(param1.points,param2.points)
	names(params.points)<-focal.param
	params.points.list<-split(params.points,row(params.points),drop=TRUE)
	likelihoods<-(sapply(params.points.list,optimizeSemifixed,phy=phy,data=data, model=model, simmap.tree=simmap.tree,scaleHeight=scaleHeight,root.station=root.station,lb=lb,ub=ub,clade=clade,mserr=mserr,globalMLE=globalMLE,simplify=TRUE))
	
	#include the MLE in the set
	likelihoods<-c(likelihoods,globalMLE$loglik)
	param1.points<-c(param1.points,as.numeric(focal.param.df[3,1]))
	param2.points<-c(param2.points,as.numeric(focal.param.df[3,2]))
	
	likelihoods.rescaled<-(-1)*likelihoods
	likelihoods.rescaled<-likelihoods.rescaled-min(likelihoods.rescaled)
	goodLikelihoodsIndex<-which(likelihoods.rescaled<likelihood.boundary)
	while(length(goodLikelihoodsIndex)<5){ #just to deal with extreme cases with few points
		likelihood.boundary<-2*likelihood.boundary
		goodLikelihoodsIndex<-which(likelihoods.rescaled<likelihood.boundary)
	}
	#interpolated.points<-interp(x= param1.points[which(likelihoods.rescaled<likelihood.boundary)],y=param2.points[which(likelihoods.rescaled<likelihood.boundary)],z=likelihoods.rescaled[which(likelihoods.rescaled<likelihood.boundary)],xo=seq(min(param1.points), max(param1.points), length = 400),yo=seq(min(param1.points), max(param1.points), length = 400),duplicate="median",linear=FALSE,extrap=TRUE)
	interpolated.points<-interp(x= param1.points[goodLikelihoodsIndex],y=param2.points[goodLikelihoodsIndex],z=likelihoods.rescaled[goodLikelihoodsIndex],linear=FALSE,extrap=TRUE,xo=seq(min(xlim), max(xlim), length = 400),yo=seq(min(ylim), max(ylim), length = 400))
	#image(interpolated.points$x,interpolated.points$y,log(abs(interpolated.points$z)),xlab=focal.param[1],ylab=focal.param[2],xlim=xlim,ylim=ylim,col=gray((701:1000)/1000)) #the log here is just to flatten out the space to make plotting prettier
	contour(interpolated.points, xlim=xlim,ylim=ylim,xlab=focal.param[1],ylab=focal.param[2], levels=levels,add=FALSE,lwd=lwd,...)
	if(strsplit(focal.param,"_")[[1]][1] == strsplit(focal.param,"_")[[2]][1]) {
		plot.range<-range(c(param1.points,param2.points))
		lines(x= plot.range,y= plot.range,lty="dotted")
	}
	points(x=as.numeric(focal.param.df[3,1]),y=as.numeric(focal.param.df[3,2]),pch=20,col="red")
	lines(x=c(as.numeric(focal.param.df[3,1])-1.96*as.numeric(focal.param.df[4,1]),as.numeric(focal.param.df[3,1])+1.96*as.numeric(focal.param.df[4,1])),y=rep(as.numeric(focal.param.df[3,2]),2))
	lines(x=rep(as.numeric(focal.param.df[3,1]),2),y=c(as.numeric(focal.param.df[3,2])-1.96*as.numeric(focal.param.df[4,2]),as.numeric(focal.param.df[3,2])+1.96*as.numeric(focal.param.df[4,2])))
	
	finalResults<-data.frame(param1.points,param2.points,likelihoods)
	names(finalResults)<-c(focal.param,"loglik")
	return(finalResults)             
}
