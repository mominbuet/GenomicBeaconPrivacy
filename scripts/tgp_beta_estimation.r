source("~/suyashs_newprojects/beacon/scripts/utils.r")

args<-commandArgs(T)
library("data.table")
outfile<-args[[1]]
fnamelist=args[2:length(args)]
ans=c()
for(fname in fnamelist){
	ob=data.frame(fread(fname))
	pop=strsplit(basename(fname),'[.]')[[1]][1]
	betaparams=round(estimatebeta(ob[,dim(ob)[2]]),digits=4)
	ans=rbind(ans,c(pop,betaparams))
}

write.table(ans,file=outfile,quote=F,row.names=F,col.names=F,sep=",")
