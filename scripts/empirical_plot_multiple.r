args<-commandArgs(T)

#Plot results for multiple analysis (CEU, CEU+JPT, CEU+YRI)

source("~/suyashs_newprojects/beacon/scripts/utils.r")

print(args)

outfile<-args[[1]]

filelist<-args[seq(2,length(args),2)]
namelist<-args[seq(3,length(args),2)]
namelist<-gsub(" ","\n",namelist)

cat(filelist,"\n")
cat(namelist,"\n")

df=data.frame(do.call(rbind,mapply(function(f,n){x=read.table(f,header=T);nrows=dim(x)[1];x$BeaconPop=rep(n,nrows);x },filelist,namelist,SIMPLIFY=F)))

head(df)

library("ggplot2")
library("reshape2")
library("scales")
library("grid")

#theorycopy=df
#theorycopy$Power=apply(theorycopy,1,function(x){theorypower(n=x[1],N=65,delta=1e-10)})
#df$Type<-rep("Empirical",dim(df)[1])
#theorycopy$Type<-rep("Theory",dim(theorycopy)[1])

#combined=rbind(df,theorycopy)

p<-ggplot(df,aes(x=Nsnps,y=Power,colour=BeaconPop))+geom_line(lwd=1.5)+ scale_x_log10()+xlab("Number of SNPs")
p<-p+scale_colour_manual(values=cbbPalette,name="Beacon\nDataset")
#p<-ggplot(combined,aes(x=Nsnps,y=Power,linetype=Type))+geom_line(lwd=1.5)+ scale_x_log10()+xlab("Number of SNPs")
png(outfile,width=600)
print(enhance(p))
dev.off()
svg(gsub(".png",".svg",outfile),width=10,height=9)
print(enhance(p))
dev.off()
