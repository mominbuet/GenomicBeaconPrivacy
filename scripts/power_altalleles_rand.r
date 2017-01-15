args<-commandArgs(T)

source("/import/helium-share/csgrad/azizmma/bustamante/scripts/utils.r")

inprefix<-args[[1]]
inno<-as.numeric(args[[2]])
outprefix<-args[[3]]
outno<-as.numeric(args[[4]])
nsnps<-as.numeric(args[[5]])
fpvec<-as.numeric(args[6:length(args)])

nreps=30

inmat=counthitsrand(inprefix,inno,nreps,nsnps)
outmat=counthitsrand(outprefix,outno,nreps,nsnps)

for(fp in fpvec){
power=binomialpower2(inmat,outmat,fp)
#print(power)
cat(nsnps,fp,power,"\n")
}
