args<-commandArgs(T)

related_alleles<-function(g,f,phi){
	alleles=c(floor(g/2),ceiling(g/2))
	sum(sapply(alleles,function(x){ifelse(runif(1)<=phi,x,rbinom(1,size=1,prob=f))}))
}

make_relative<-function(x,phi,fvec){
	mapply(FUN=related_alleles,x,fvec,phi)
}

generate_relateds<-function(df,relatedness,freqs){
	t(apply(df,1,function(x){make_relative(x,relatedness,freqs)}))
}

#set.seed(2)

frqfile<-args[[1]]

relatedness<-args[[2]]

batchfiles<-args[3:length(args)]

print(batchfiles)

library("data.table")

freqs=scan(frqfile)

for(file in batchfiles){
	
	cat(file,"\n")
	df=data.frame(fread(file))
	gc()
	relatives=generate_relateds(df,relatedness,freqs)
	outfile=paste0(file,".related_",relatedness)
	write.table(relatives,outfile,quote=F,row.names=F,col.names=F)
}
