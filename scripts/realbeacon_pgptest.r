args<-commandArgs(T)

source("~/suyashs_newprojects/beacon/scripts/utils.r")

beaconfile<-args[[1]]
pgpfile<-args[[2]]
outfile<-args[[3]]

beacondata<-read.table(beaconfile,sep=',',header=F,stringsAsFactors=F)
colnames(beacondata)<-c("Beacon","Ninds")

pgphits=read.table(pgpfile,sep=",",stringsAsFactor=F,header=F)
colnames(pgphits)<-c("Beacon","Nhits")

library("plyr")

mergedata=join(pgphits,beacondata,type="inner")
mergedata$Pval=apply(mergedata[,c("Ninds","Nhits")],1,function(x){print(x);getpval(N=x[1],n=1000,nhits=x[2])})

results=mergedata[,c("Beacon","Ninds","Nhits","Pval")]

head(results)

write.table(results,file=outfile,sep=',',quote=F,col.names=T,row.names=F)
