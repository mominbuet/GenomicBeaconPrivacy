args<-commandArgs(T)

source("~/suyashs_newprojects/beacon/scripts/utils.r")
library("ggplot2")
library("scales") 
library("grid") 

paramslist=vector("list", length(args))
idx=0
for(file in args){
	idx=idx+1
	f=read.csv(file,header=F)
	dataset=toupper(strsplit(file,split="_")[[1]][1])
	paramslist[[idx]]=c(dataset,f[1,2],f[1,3])
}

head(paramslist)

querysizes=c(1000,2000,5000,10000,20000,40000,100000)
querysizes=10^(2:6 )
beaconsizes=c(100,1000,5000,10000)
sfseffect=expand.grid(querysizes,paramslist,beaconsizes)
colnames(sfseffect)=c("Nsnps","Params","N")
sfseffect$Dataset=apply(sfseffect,1,function(x){x[2][[1]][1]})
sfseffect$alphasfs=apply(sfseffect,1,function(x){as.numeric(x[2][[1]][2])})
sfseffect$betasfs=apply(sfseffect,1,function(x){as.numeric(x[2][[1]][3])})
sfseffect$Power=apply(sfseffect,1,function(x){theorypower(n=x$Nsnps[1],N=x$N[1],alphasfs=1+x$alphasfs[1],betasfs=1+x$betasfs[1],alpha=0.05)})
p<-ggplot(sfseffect,aes(x=Nsnps,y=Power,colour=Dataset))+geom_line(lwd=1.5,position=position_jitter(w=0.15,h=0))+facet_wrap(~N,ncol=2)+theme(strip.text.x = element_text(size=20))
#p<-p+ scale_x_log10(breaks=querysizes,labels=as.character(querysizes))
p<-p+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))
p<-p+scale_colour_manual(values=cbbPalette,name="Population")+xlab("Number of SNPs")
png("sfs_effect_power.png",width = 960, height = 960)
print(enhance(p))
dev.off()
