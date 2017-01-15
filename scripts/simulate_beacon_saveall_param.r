args<-commandArgs(T)

source("/import/helium-share/csgrad/azizmma/bustamante/scripts/utils.r")

ninds<-3e3
nsnps<-as.numeric(args[[2]])
npool=as.numeric(args[[3]])
print(npool)
set.seed(1)

savedir<-args[[1]]

fvec=generatefs(npool,nsnps)
print(quantile(fvec,c(1:10)/10))
betaparams=estimatebeta(fvec)

print(betaparams)

nbatch=100
nreps=ninds/nbatch
#0.1284, 1.13096
write(fvec,file=paste0(savedir,"/Frequencies.txt"),ncol=1)
cat("Wrote frequencies ",nreps,"\n")

for(i in 1:nreps){
	cat("Batch",i,"\n")
	siminds=do.call(rbind,lapply(1:nbatch ,function(x){sapply(fvec,rbinom,n=1,size=2)} ))
	fname=paste0(savedir,"/Batch",i,".txt")
	write.table(siminds,file=fname,quote=F,row.names=F,col.names=F)
	
}
cat("Wrote inds\n")
