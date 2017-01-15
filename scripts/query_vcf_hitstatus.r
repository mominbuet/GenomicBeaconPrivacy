library("data.table")
library("plyr")

args<-commandArgs(T)

queryind=args[[1]] #"hmp3_ceu_altalleles/NA06984.altalleles"
beacondir=args[[2]] #"../../data/2014_11_24/phase1_ceu_beacon"
outfile=args[[3]]

query=data.frame(fread(queryind,stringsAsFactors=FALSE),stringsAsFactors=FALSE)
chr_cols<-c("Chr","Pos","Allele")
colnames(query)<-chr_cols

hitfreq=c()
for(chr in 1:22 ){
	cat("Chromosome",chr,"\n")
	chrfile=paste(beacondir,"/",chr,".altfreqs",sep="")
	chrtable=read.table(chrfile,stringsAsFactors=FALSE)
	colnames(chrtable)<-c(chr_cols,"AltF")
	#print(head(chrtable))
	joins=join(query[query$Chr==chr,],chrtable,by=chr_cols,type="left")
	#print(dim(joins))
	hitfreq=c(hitfreq,joins$AltF)
}
write(x=hitfreq,file=outfile,ncolumns=1)
