source("/import/helium-share/csgrad/azizmma/bustamante/scripts/utils.r")
library("data.table")

args<-commandArgs(T)

frqfile=args[[1]]
freqs=scan(frqfile)

#set.seed(3)

delta=as.numeric(args[[2]])

infiles=args[3:length(args)]

cat("Reading files now\n")
for(file in infiles){
	print(file)
	df=as.matrix(fread(file))
	cat("Distorting\n")
	#print(df)
	queryinds=distort(df,delta)
	outstr=paste0("Query_",delta,"_Batch")
	outfile=gsub("Batch",outstr,file)
	print(outfile)
	write.table(queryinds,file=outfile,quote=F,row.names=F,col.names=F)
}
