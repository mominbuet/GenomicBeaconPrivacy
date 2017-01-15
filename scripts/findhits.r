args<-commandArgs(T)

#Returns for each query position whether it was found in the beacon list of positions or not

beaconfile<-args[[1]]

beacontable=read.table(beaconfile)

beaconpos=beacontable[,1]

infiles<-args[2:length(args)]

for(file in infiles){
	print(file)
	outfile=gsub("altalleles","queryresults",file)
	qpos=scan(file)
	qres=0+(qpos %in% beaconpos)
	results=data.frame(qpos=qpos,qres=qres)
	write.table(results,file=outfile,quote=F,row.names=F,col.names=F)
}
