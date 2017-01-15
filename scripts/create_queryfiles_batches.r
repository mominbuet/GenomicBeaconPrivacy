args<-commandArgs(T)
#Create query file with extension .altalleles from the genotype matrix

outprefix<-args[[1]]

startidx<-as.numeric(args[[2]])

infiles<-args[3:length(args)]

library("data.table")

idx=startidx
for(file in infiles){
	df=as.matrix(fread(file))
	ninds=dim(df)[1]
	cat("ninds ",ninds)
	for(i in 1:ninds ){
		outfile=paste0(outprefix,idx,".altalleles")
		print(outfile)
		idx=idx+1	
		querypos=which(df[i,]==1)
		write(querypos,file=outfile,ncolumn=1)
	}
}
