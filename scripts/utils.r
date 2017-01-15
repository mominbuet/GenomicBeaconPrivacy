# Colorblind-friendly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Flips proportion delta of heterozygous SNPs per individual (rows of xmat), and replaces the same number of hom.refs by hets to create a noisy copy
distort<-function(xmat,delta){
	xd=xmat
	for(i in 1:dim(xmat)[1]){
				orig=xmat[i,]
				flipped=sapply(orig,function(x){ifelse(runif(1)<=delta && x==1,0,x)})
				nflipped=sum(orig!=flipped)
				flipped[sample(which(orig==0))[1:nflipped]]=1
				xd[i,]=flipped
			}
	xd
}

#Method of moments estimator of beta
estimatebeta<-function(hitf){
	meanf=mean(hitf)
	varf=var(hitf)

	alphaest=meanf*(meanf*(1-meanf)/varf-1)
	betaest=(1-meanf)*(meanf*(1-meanf)/varf-1)
	c(alphaest,betaest)
}

#Generate frequencies according to the standard neutral model for a given population size
generatefs<-function(npool,nsnps){
	countseq=1:(npool-1)
	sample(countseq,nsnps,replace=T,prob=1/countseq)/npool
}

#Compute D_{N-1} as described in the paper
dvalue<-function(N,alphasfs,betasfs){
  Dval=prod(sapply(0:(2*N-3),function(x){(betasfs+x)/(alphasfs+betasfs+x)}))
  Dval
}

bvalue<-function(N,alphasfs,betasfs,delta){
	Dval=dvalue(N,alphasfs,betasfs)
	Dp=dvalue(N+1,alphasfs,betasfs)
	B=log(Dp)-log(delta*Dval)
	B
}

cvalue<-function(N,alphasfs,betasfs,delta){
	D=dvalue(N,alphasfs,betasfs)
	Dp=dvalue(N+1,alphasfs,betasfs)
	C=log(delta*D*(1-Dp))-log(Dp*(1-delta*D))
	C
}

#power according to the normal approximation
llrpower<-function(alpha=0.05,n,N,delta=0.001,alphasfs=1.01,betasfs=2.01){
	Dn=dvalue(N,alphasfs,betasfs)
	Dnp=dvalue(N+1,alphasfs,betasfs)
	zalpha=qnorm(alpha)
	zbeta=(-sqrt(n)*(delta*Dn-Dnp)+zalpha*sqrt(Dnp*(1-Dnp)))/(sqrt(delta*Dn*(1-delta*Dn)))
	pnorm(zbeta)
}

#Power according to theory using the Binomial expressions
theorypower<-function(alpha=0.05,n,N,delta=0.001,alphasfs=1.01,betasfs=2.01){
	#cat("Alpha is",alpha,"\n")
	Dn=dvalue(N,alphasfs,betasfs)
	Dnp=dvalue(N+1,alphasfs,betasfs)
	thresh=qbinom(alpha,size=n,prob=1-Dnp,lower.tail=F)
	pbinom(thresh,size=n,prob=1-delta*Dn,lower.tail=F)	
}

relativetheorypower2<-function(alpha=0.05,n,N,delta=0.001,alphasfs=1.01,betasfs=2.01,phi=1){
	Dn=dvalue(N,alphasfs,betasfs)
	Dnp=dvalue(N+1,alphasfs,betasfs)
	Dnint=dvalue(N+0.5,alphasfs,betasfs)
	#cat("Vals are",Dn,Dnp,Dnint,"\n")
	nullprob=Dnp
	altprob=delta*Dn+(1-phi)^2*(1-2*delta)*Dnp+1*phi*(1-phi)*(1-2*delta)*Dnint
	thresh=qbinom(alpha,size=n,prob=nullprob,lower.tail=T)
	#cat("n and thresh and nullprob and altprob are",n,thresh,nullprob,altprob,"\n")
	pbinom(thresh,size=n,prob=altprob,lower.tail=T)	
}

#Power for detecting relatives
relativetheorypower<-function(alpha=0.05,n,N,delta=0.001,alphasfs=1.01,betasfs=2.01,phi=1){
	Dn=dvalue(N,alphasfs,betasfs)
	Dnp=dvalue(N+1,alphasfs,betasfs)
	Dnint=dvalue(N+0.5,alphasfs,betasfs)
	#cat("Vals are",Dn,Dnp,Dnint,"\n")
	nullprob=1-Dnp
	altprob=1-delta*Dn-(1-phi)^2*(1-2*delta)*Dnp-1*phi*(1-phi)*(1-2*delta)*Dnint
	thresh=qbinom(alpha,size=n,prob=nullprob,lower.tail=F)
	#cat("n and thresh and nullprob and altprob are",n,thresh,nullprob,altprob,"\n")
	pbinom(thresh,size=n,prob=altprob,lower.tail=F)	
}

theorypower2<-function(alpha=0.05,n,N,delta=0.001,alphasfs=0.12,betasfs=1.12){
	Dn=2*alphasfs*exp(lgamma(betasfs+2*N-1)+lgamma(alphasfs+betasfs)-lgamma(betasfs)-lgamma(alphasfs+betasfs+2*N)) #dvalue(N,alphasfs,betasfs)
	Dnp=2*alphasfs*exp(lgamma(betasfs+2*N+1)+lgamma(alphasfs+betasfs)-lgamma(betasfs)-lgamma(alphasfs+betasfs+2*N+2))#dvalue(N+1,alphasfs,betasfs)
	oldDn=dvalue(N,alphasfs,betasfs)
	oldDnp=dvalue(N+1,alphasfs,betasfs)
	cat("old and new",Dn,Dnp,oldDn,oldDnp,"\n")
	thresh=qbinom(alpha,size=n,prob=1-Dnp,lower.tail=F)
	pbinom(thresh,size=n,prob=1-delta*Dn,lower.tail=F)	
}

getpower<-function(indata,alpha){
	#cat("Alpha is",alpha,"\n")
	fpthresh=quantile(indata[indata$In==0,2],alpha)
	power=sum(indata[indata$In==1,2]<fpthresh)/sum(indata$In==1)
}

#This has a bad transposition, reps are in rows,not columns. each column is one ind. so thresholds need to be found across rows, not columns
binomialpower<-function(inmat,outmat,alpha){
	fpthresh=apply(outmat,2,function(x){quantile(x,1-alpha)})
	pvec=rep(0,dim(inmat)[2])
	for(i in 1:length(pvec)){
		pvec[i]=sum(inmat[,i]>fpthresh[i])/dim(inmat)[1]
	}
	mean(pvec)
}

#Power calculation for many replicates of the test
#inmat=Test statistic for individuals known to be in beacon
#outmat=Test statistic for individuals known to be outside beacon
#Both have replicates in rows, individuals in columns
binomialpower2<-function(inmat,outmat,alpha){
	#cat(alpha)
	fpthresh=apply(outmat,1,function(x){quantile(x,1-alpha)})
	pvec=rep(0,dim(inmat)[1])
	for(i in 1:length(pvec)){
		pvec[i]=sum(inmat[i,]>fpthresh[i])/dim(inmat)[2]
	}
	mean(pvec)
}

#Given the number of hits for a whole genome, subsample to the number of hits for a smaller set of queries.
#inprefix=File path ; inno=Individual number ; nreps=No. of replicates ; nsnps=No. of queries
counthits<-function(inprefix,inno,nreps,nsnps){
	inmat=matrix(0,nrow=nreps,ncol=inno)
	for(i in 1:inno ){
		fname=paste0(inprefix,i,".queryresults")
		status=read.table(fname)[,2]
		inmat[,i]=replicate(nreps,simplify=T,expr=sum(sample(status,size=nsnps)))
	}
	#print(inmat)
	inmat
}
counthitsrand<-function(inprefix,inno,nreps,nsnps){
	inmat=matrix(0,nrow=nreps,ncol=inno)
	for(i in 1:inno ){
		fname=paste0(inprefix,i,".rand_queryresults")
		status=read.table(fname)[,2]
		inmat[,i]=replicate(nreps,simplify=T,expr=sum(sample(status,size=nsnps)))
	}
	#print(inmat)
	inmat
}
#Better formatting of ggplot axis labels with a 45-degree angle of labels
enhance45<-function(p){
p<-p+theme_bw()
p<-p+theme(axis.title.x = element_text(face="bold", size=20),
           axis.text.x  = element_text(size=16,angle=45,vjust=0.5,colour='black'))
p<-p+theme(axis.title.y = element_text(face="bold", size=20),
           axis.text.y  = element_text(size=16,colour='black'))
p<-p+ theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 16))
p
}

#Better formatting of ggplot axis labels 
enhance<-function(p,xangle=0){
p<-p+theme_bw()
p<-p+theme(axis.title.x = element_text(face="bold", size=20),
           axis.text.x  = element_text(size=16,angle=xangle,vjust=0.5,colour='black'))
p<-p+theme(axis.title.y = element_text(face="bold", size=20),
           axis.text.y  = element_text(size=16,colour='black'))
p<-p+ theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 16))
p<-p+theme(panel.margin = unit(1, "cm")) 
p
}

#Beacon size required to ensure limited power against n queries. Uses approximation
approxN<-function(n,alpha=0.05,beta=0.95,delta=0.001,alphasfs=1.2,betasfs=2.1){
	zalpha=qnorm(alpha)
	zbeta=qnorm(beta)
	ans=((n*gamma(alphasfs+betasfs)/(gamma(betasfs)*(zalpha-zbeta*delta)^2))^(1/alphasfs))/2
}

#Beacon size required to ensure limited power against n queries.
exactN<-function(n,alpha=0.05,beta=0.95,delta=0.001,alphasfs=1.2,betasfs=2.1){
	zalpha=qnorm(alpha)
	zbeta=qnorm(beta)
	frac=1-((n*delta-zbeta^2*delta-2*zbeta*zalpha*sqrt(delta))/(n+zalpha^2))
	ans=alphasfs/frac-alphasfs-betasfs+1
	ans=ans/2	
}

#No of queries required to identify a relative (relatedness=phi) in a beacon of size N
minn_phi<-function(N,alpha=0.05,beta=0.95,delta=1e-2,alphasfs=1.0,betasfs=2.0,phi=1){
	Dn=dvalue(N,alphasfs,betasfs)
	Dnp=dvalue(N+1,alphasfs,betasfs)
	Dnint=dvalue(N+0.5,alphasfs,betasfs)
	#cat("Vals are",Dn,Dnp,Dnint,"\n")
	nullprob0=Dnp
	altprob0=delta*Dn+(1-phi)^2*(1-2*delta)*Dnp+1*phi*(1-phi)*(1-2*delta)*Dnint
	zalpha=qnorm(alpha)
	zbeta=qnorm(beta)
	n=((zalpha*sqrt(nullprob0*(1-nullprob0))-zbeta*sqrt(altprob0*(1-altprob0)))/(altprob0-nullprob0))^2	
	ceiling(n)
}

#Theoretical p-value calculation for the PGP experiment
getpval<-function(N,n,nhits,alphasfs=1.0,betasfs=2.0){
	Dnp=dvalue(N+1,alphasfs,betasfs)
	pbinom(n-nhits,size=n,prob=Dnp) # Count number of misses instead of hits
}

#Power for the case when some variation is censored. Not used.
censoredpower<-function(alpha=0.05,n,N,delta=0.001,alphasfs=1.01,betasfs=2.01,fthresh){
	Dn=dvalue(N,alphasfs,betasfs)
	Dnp=dvalue(N+1,alphasfs,betasfs)
	Enminusone=pbeta(fthresh,alphasfs,betasfs)+delta*Dn*(1-pbeta(fthresh,alphasfs,betasfs+2*N-2))
	Fn=pbeta(fthresh,alphasfs,betasfs)+Dnp*(1-pbeta(fthresh,alphasfs,betasfs+2*N))
	zalpha=qnorm(alpha)
	zbeta=(-sqrt(n)*(Enminusone-Fn)+zalpha*sqrt(Fn*(1-Fn)))/(sqrt(Enminusone*(1-Enminusone)))
  #  zbeta=((1-delta)*sqrt(n*D)+zalpha*sqrt(1-D))/sqrt(delta*(1-delta*D))
  #cat("zbeta is",zbeta,"\n")
	pnorm(zbeta)
}
