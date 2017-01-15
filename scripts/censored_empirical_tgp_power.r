args<-commandArgs(T)

#Script for generating results from censoring, modifies how singleton responses are processed

source("~/suyashs_newprojects/beacon/scripts/utils.r")

createmat<-function(idlist,nreps,nsnps){
	ninds=length(idlist)
	mat=matrix(0,nrow=nreps,ncol=ninds)
	for(ino in 1:ninds ){
		#cat(ino,"\n")
		hitlist=scan(paste0("./tgp_ceu_query_altalleles/",idlist[ino],".hits"),quiet=T)
		hitlist1=ifelse(is.na(hitlist) | hitlist<= 1/(2*ninds),0,1) #Added extra term for censoring output for singletons
		mat[,ino]=replicate(nreps,simplify=T,expr=sum(sample(hitlist1,size=nsnps)))
	}
	mat
}

beaconstatusfile<-args[[1]]
nsnps<-as.numeric(args[[2]])

beaconstatus=read.table(beaconstatusfile,stringsAsFactors=F)


nreps=30
ids_inbeacon=beaconstatus[beaconstatus[,2]==1,1]
ids_outbeacon=beaconstatus[beaconstatus[,2]==0,1]

#cat("inmat\n")
inmat=createmat(ids_inbeacon,nreps,nsnps)
#cat("outmat\n")
outmat=createmat(ids_outbeacon,nreps,nsnps)

alpha=0.05
power=binomialpower2(inmat,outmat,alpha)

cat(nsnps,alpha,power,"\n")
