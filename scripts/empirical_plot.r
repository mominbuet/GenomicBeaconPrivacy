args<-commandArgs(T)
#Plot results for a single empirical analysis

source("~/suyashs_newprojects/beacon/scripts/utils.r")

infile<-args[[1]]

outfile<-args[[2]]

df=read.table(infile,header=T)

library("ggplot2")
library("reshape2")
library("scales")
library("grid")

theorycopy=df
theorycopy$Power=apply(theorycopy,1,function(x){theorypower(n=x[1],N=65,delta=1e-10)})
df$Type<-rep("Empirical",dim(df)[1])
theorycopy$Type<-rep("Theory",dim(theorycopy)[1])

combined=rbind(df,theorycopy)

p<-ggplot(combined,aes(x=Nsnps,y=Power,linetype=Type))+geom_line(lwd=1.5)+ scale_x_log10()+xlab("Number of SNPs")
png(outfile,width=600)
print(enhance(p))
dev.off()
