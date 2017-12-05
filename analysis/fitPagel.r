## function fits Pagel '94 model of correlated evolution of two binary characters
## uses ape::ace or geiger::fitDiscrete internally
## written by Liam J. Revell 2014



library(shape)


fitPagel<-function(tree,x,y,...){
	if(hasArg(method)) method<-list(...)$method
	else method<-"ace"
	if(method=="fitDiscrete"){
		#chk<-.check.pkg("geiger")
		chk = T
		if(!chk){
			cat("  method = \"fitDiscrete\" requires the package \"geiger\"\n")
			cat("  Defaulting to method = \"ace\"\n\n")
			method<-"ace"
			fitDiscrete<-function(...) NULL
		}
	}
	if(!is.factor(x)) x<-as.factor(x)
	levels.x<-levels(x)
	if(!is.factor(y)) y<-as.factor(y)
	levels.y<-levels(y)
	y<-y[names(x)]
	if(length(levels.x)!=2||length(levels.y)!=2)
		stop("Only binary characters for x & y currently permitted.")
	xy<-setNames(factor(paste(x,y,sep="|"),
		levels=sapply(levels.x,paste,levels.y,sep="|")),
		names(x))
	## fit independent model

	## model shapes
	if(!hasArg(iQ)){
	iQ<-matrix(c(0,1,2,0,3,0,0,2,4,0,0,1,0,4,3,0),4,4,byrow=TRUE)
	}
	rownames(iQ)<-colnames(iQ)<-levels(xy)
	if(!hasArg(dQ)){
	dQ<-matrix(c(0,1,2,0,3,0,0,4,5,0,0,6,0,7,8,0),4,4,byrow=TRUE)
	}
	rownames(dQ)<-colnames(dQ)<-levels(xy)
	
	# base degrees of freedom on difference in number of rates
	# (original code does it based on number of levels)
	dfx=length(levels(as.factor(dQ[dQ!=0]))) - length(levels(as.factor(iQ[iQ!=0])))
	

	fit.iQ<-if(method=="fitDiscrete") fitDiscrete(tree,xy,model=iQ) else ace(xy,tree,type="discrete",model=iQ)
	## fit dependendent model
	fit.dQ<-if(method=="fitDiscrete") fitDiscrete(tree,xy,model=dQ) else ace(xy,tree,type="discrete",model=dQ)
	## back translate independent model
	if(method=="fitDiscrete") iQ<-geiger:::.Qmatrix.from.gfit(fit.iQ)
	else {
		I<-fit.iQ$index.matrix
		I[I==0]<-NA
		iQ<-apply(I,2,function(i,x) x[i],x=fit.iQ$rates)
		iQ[is.na(iQ)]<-0
		diag(iQ)<--rowSums(iQ)
		rownames(iQ)<-colnames(iQ)
	}
	## dependent model
	if(method=="fitDiscrete") dQ<-geiger:::.Qmatrix.from.gfit(fit.dQ)
	else {
		I<-fit.dQ$index.matrix
		I[I==0]<-NA
		dQ<-apply(I,2,function(i,x) x[i],x=fit.dQ$rates)
		dQ[is.na(dQ)]<-0
		diag(dQ)<--rowSums(dQ)
		rownames(dQ)<-colnames(dQ)
	}
	## assemble object to return
	obj<-list(independent.Q=iQ,
		dependent.Q=dQ,
		independent.logL=logLik(fit.iQ),
		dependent.logL=logLik(fit.dQ),
		lik.ratio=2*(logLik(fit.dQ)-logLik(fit.iQ)),
		P=pchisq(2*(logLik(fit.dQ)-logLik(fit.iQ)),
		#df=length(levels(x))+length(levels(y)),
		df=dfx,
		lower.tail=FALSE))
	class(obj)<-"fitPagel"
	obj
}

print.fitPagel<-function(x,...){
	cat("\n  Pagel's binary character correlation test:\n")
	cat("\nIndepedent model rate matrix:\n")
	print(x$independent.Q)
	cat("\nDependent model rate matrix:\n")
	print(x$dependent.Q)
	cat("\nModel fit:\n")
	obj<-matrix(c(x$independent.logL,x$dependent.logL),2,1)
	rownames(obj)<-c("independent","dependent")
	colnames(obj)<-"log-likelihood"
	print(obj)
	cat("\nHypothesis test result:\n")
	cat(paste("  likelihood-ratio: ",signif(x$lik.ratio,7),"\n"))
	cat(paste("  p-value: ",signif(x$P,7),"\n"))
}


DrawRateGraph <-function(fit.Discrete){
	
	minArrowWidth = 1
	maxArrowWidth = 10
	
	q = fit.Discrete$dependent.Q
	
	range = q[upper.tri(q)|lower.tri(q)]
	q = (q-min(range)) /max(range-min(range))
	
	plot(-1:5,-1:5,xaxt='n',yaxt='n',bty='n',col='white',xlab='',ylab='')
	dx = 1.5
	rNames = rownames(q)
	stateNames = gsub("\\|","\n",rownames(q))
	pos = list(c(1,1),c(1,3),c(3,1),c(3,3))
	for(p in 1:length(pos)){
		text(pos[p][[1]][1],pos[p][[1]][2],stateNames[p])
	}
	
	for(i in 1:length(rNames)){
		for(j in 1:length(rNames)){
			if(i!=j){
				from = rNames[i]
				to = rNames[j]
				strength = minArrowWidth + (q[from,to]*maxArrowWidth)
				loc1 = pos[i][[1]]
				loc2 = pos[j][[1]]
				
				if((loc1[1]==loc2[1]) | (loc1[2]==loc2[2])){
				
				

				if(loc1[1]==loc2[1]){
					adj = (loc1[1] - 2)
					if(i>j){
						loc1[1] = loc1[1]+adj
						loc2[1] = loc2[1]+adj
					} else{
						loc1[1] = loc1[1]+(adj*dx)
						loc2[1] = loc2[1]+(adj*dx)
					}
				} 
				if(loc1[2]==loc2[2]){
					adj = (loc1[2] - 2)
					if(i>j){
						loc1[2] = loc1[2]+adj
						loc2[2] = loc2[2]+adj
					} else{
						loc1[2] = loc1[2]+(adj*dx)
						loc2[2] = loc2[2]+(adj*dx)
					}
				}
				
				
				arrows(loc1[1],loc1[2],loc2[1],loc2[2],lwd=strength)
				
			}
			}
		}
	}
	
}



DrawRateGraph2 <-function(fit.Discrete,arrow.cols=rep(1,8),independent=F){
	
	minArrowWidth = 1
	maxArrowWidth = 10
	q = fit.Discrete$dependent.Q
	if(independent){
		q = fit.Discrete$independent.Q
	}
	
	range = q[upper.tri(q)|lower.tri(q)]
	if(independent){
		q2 = fit.Discrete$dependent.Q
		range = q2[upper.tri(q2)|lower.tri(q2)]
	}
	q = (q-min(range)) /max(range-min(range))
	
	par(mar=c(0,0,0,0))
	plot(-2:6,-2:6,xaxt='n',yaxt='n',bty='n',col='white',xlab='',ylab='')
	
	dx = 1.5
	rNames = rownames(q)
	stateNames = gsub("\\|","\n",rownames(q))
	pos = list(c(1,1),c(1,3),c(3,1),c(3,3))
	pos.text = list(c(0,0),c(0,4),c(4,0),c(4,4))
	for(p in 1:length(pos)){
		text(pos.text[p][[1]][1],pos.text[p][[1]][2],stateNames[p])
	}
	cx = 1
	for(i in 1:length(rNames)){
		for(j in 1:length(rNames)){
			if(i!=j){
				from = rNames[i]
				to = rNames[j]
				strength = minArrowWidth + (q[from,to]*maxArrowWidth)
				loc1 = pos[i][[1]]
				loc2 = pos[j][[1]]
				
				if((loc1[1]==loc2[1]) | (loc1[2]==loc2[2])){
				
				

				if(loc1[1]==loc2[1]){
					adj = (loc1[1] - 2)
					if(i>j){
						loc1[1] = loc1[1]+adj
						loc2[1] = loc2[1]+adj
					} else{
						loc1[1] = loc1[1]+(adj*dx)
						loc2[1] = loc2[1]+(adj*dx)
					}
				} 
				if(loc1[2]==loc2[2]){
					adj = (loc1[2] - 2)
					if(i>j){
						loc1[2] = loc1[2]+adj
						loc2[2] = loc2[2]+adj
					} else{
						loc1[2] = loc1[2]+(adj*dx)
						loc2[2] = loc2[2]+(adj*dx)
					}
				}
				
				
	par(ljoin=3)
	Arrows(loc1[1],loc1[2],loc2[1],loc2[2],lwd=strength,col=arrow.cols[cx], arr.type='triangle')
	par(ljoin=0)
				cx = cx + 1
			}
			}
		}
	}
	
}



DrawRateGraphWithHist = function(s,fitDiscrete,drawAllHist=F,arrow.cols=c(1,5,2,4,3,1,3,2),box.cols = c(2,3,1,5,4),independent=F, n.breaks=10){
	pos.h = list(c(),c(0.2,2),c(1.5,1),c(),c(-1.8,2),c(),c(),c(),c(1.5,-1.2),c(),c(),c(),c(),c(1.5,5.2),c(),c())
	
	if(drawAllHist){
		pos.h = list(c(),c(0.2,2),c(1.5,1),c(),c(-1.8,2),c(),c(),c(1.5,3),c(1.5,-1.2),c(),c(),c(2.8,2),c(),c(1.5,5.2),c(4.8,2),c())
	
	}
	
	breaks = c(seq(0,1,length.out=n.breaks))
	if(independent & sum(arrow.cols==c(1,5,2,4,3,1,3,2))==8){
		arrow.cols = c(1,4,2,4,3,1,3,2)
		box.cols = c(2,3,1,4,4)
	}
	
	
	DrawRateGraph2(fitDiscrete,arrow.cols,independent)
	cx = 1
	for(i in 1:16){
			if(length(pos.h[[i]])==2){
			if(!independent | cx!=5){
			h = hist(s[,i],plot=F,breaks=breaks)
			y =h$density
			y = 0.98 * ((y+abs(min(y)))/max(y))
			y.trans= y+(pos.h[[i]][2])-0.5
			
			x = seq(pos.h[[i]][1],pos.h[[i]][1]+1,length.out=length(y))
			lines(c(pos.h[[i]][1]+0.5,pos.h[[i]][1]+0.5),c(pos.h[[i]][2]+0.5,pos.h[[i]][2]-0.5),col='gray')
			lines(x,y.trans)
			rect(pos.h[[i]][1],pos.h[[i]][2]-0.5,pos.h[[i]][1]+1,pos.h[[i]][2]+0.5,border=box.cols[cx])

		}
		cx= cx+1
		}
	}
}



