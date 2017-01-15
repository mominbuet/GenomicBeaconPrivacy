args<-commandArgs(T)

library("data.table")

outfile<-args[[1]]
print(outfile)

infiles<-args[2:length(args)]
print(infiles)

altallelecounts=Reduce("+",lapply(infiles,function(x){colSums(as.matrix(fread(x)))}))
#altallelecounts=4000000
print(length(altallelecounts))

beacon_positions=which(altallelecounts>0)
print(length(beacon_positions))
beacon_counts=altallelecounts[beacon_positions]

df=data.frame(pos=beacon_positions,cts=beacon_counts)

write.table(df,file=outfile,row.names=F,col.names=F,quote=F)
