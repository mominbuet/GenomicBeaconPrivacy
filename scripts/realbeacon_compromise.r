args<-commandArgs(T)

source("~/suyashs_newprojects/beacon/scripts/utils.r")

infile<-args[[1]]
outfile<-args[[2]]

beacondata<-read.table(infile,sep=',',header=T,stringsAsFactors=F)
ncol=dim(beacondata)[2]
nonagg=beacondata[beacondata$Aggregator=="No" & !is.na(beacondata[,ncol]),]
head(nonagg)

nonagg$Identical<-sapply(nonagg[,ncol],minn_phi)
nonagg$First<-sapply(nonagg[,ncol],minn_phi,phi=0.5)
nonagg$Second<-sapply(nonagg[,ncol],minn_phi,phi=0.25)
results=nonagg[,-which(names(nonagg) %in% c("Aggregator","Individual.data"))]

head(results)

write.table(results,file=outfile,sep=',',quote=F,col.names=T,row.names=F)
