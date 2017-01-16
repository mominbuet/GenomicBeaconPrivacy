#args<-commandArgs(T)

rm(list = ls())
setwd("E:/Privacy Code Works/beacon_distribute/beacon_distribute")
bias<-c(50,75,90)
fprs<-c(0.1,0.1,0.1,0.1)
i=1;
source("E:\\Privacy Code Works\\beacon_distribute\\beacon_distribute\\scripts\\utils.r")
library("ggplot2")
library("scales") 
library("grid") 
library(reshape2)

fixdeltaphifp=NULL
theorycopy=NULL
#infile<-args[[1]]

# infile<-paste("E:\\Privacy Code Works\\beacon_distribute\\beacon_distribute\\allresultsElim",bias[i],".dump", sep = '')
# df=read.table(infile,header=T)


snpvec=c(10^(1:6 ),3e6)
#Nvec=sapply(snpvec,approxN,delta=1e-6)
#df1=data.frame(Nsnps=snpvec,N=Nvec)


# fixdeltaphifp<-df[df$Phi==1 & df$Delta==.01 & df$FPR==0.1 & df$Nsnps<=300000 & df$Nsnps>=10000,]
# fixdeltaphifp$Bias<-rep(paste("",bias[i]),dim(fixdeltaphifp)[1])
#
# theorycopy<-fixdeltaphifp
# theorycopy$Power=apply(theorycopy,1,function(x){llrpower(n=x[1],N=1000,delta=1e-6,alpha=0.05)})


# theorycopy$Type<-rep("Theory",dim(theorycopy)[1])

for(j in 1:length(bias)) {
  cat("",bias[j])
  infile<-paste("E:\\Privacy Code Works\\beacon_distribute\\beacon_distribute\\allresultsRand",bias[j],".dump", sep = '')
  df=read.table(infile,header=T)
  tmp<-df[df$Phi==1 & df$Delta==.01 & df$FPR==fprs[j] & df$Nsnps<=300000,]
  tmp$Bias<-rep(paste("",bias[j]),dim(tmp)[1])
  
  fixdeltaphifp<-rbind( fixdeltaphifp,  tmp )
  #theorycopy<-fixdeltaphifp
  
}
combined=rbind(fixdeltaphifp)
#}


p <- ggplot(combined,aes(x=Nsnps,y=Power,shape=Bias))+ scale_x_continuous(labels = comma)+geom_line(lwd=.8)+geom_point(size=5) +
  theme_bw() +
  theme(legend.position="top",
        legend.text = element_text(size = 16),
        panel.grid.major = element_line(color="grey", size = .2),
        panel.grid.minor = element_line(color="grey", size = .05),
        panel.border = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = .5),
        axis.line.y = element_line(color="black", size = .5),
        text = element_text(size=18))+
  ylim(0,1)+  xlab("Number of SNPs") 
ggsave(filename = "E:\\Privacy Code Works\\beacon_distribute\\beacon_distribute\\simplepowerrandAll.png",plot=p,height = 8,width=12, dpi = 400)

# print(enhance(p))
# dev.off(2)





#phi analysis
infile<-paste("E:\\Privacy Code Works\\beacon_distribute\\beacon_distribute\\allresultsRand90.dump", sep = '')
df=read.table(infile,header=T)
breakvec=c(0.01,0.025,0.05,0.1)
labelvec=c("0.01","0.025","0.05","0.1")

fixdeltansnpsphi=df[df$Nsnps==125000 & df$Delta==1e-6 & df$Phi==1 & df$FPR<=0.1, ]
theorycopy1=fixdeltansnpsphi
theorycopy1$Power=apply(theorycopy1,1,function(x){theorypower(n=5000,N=1000,delta=1e-6,alpha=x[2])})
fixdeltansnpsphi$Type<-rep("Bias 90",dim(fixdeltansnpsphi)[1])
theorycopy1$Type<-rep("Theory",dim(theorycopy1)[1])

combined=rbind(fixdeltansnpsphi,theorycopy1)

p<-ggplot(combined,aes(x=FPR,y=Power,linetype=Type))+geom_line(lwd=1)+ylim(0,1) + scale_x_log10(breaks=breakvec,labels=labelvec) +xlab("False Positive Rate") 

png(paste("E:\\Privacy Code Works\\beacon_distribute\\beacon_distribute\\fp_vs_power90.png",sep=''))
print(enhance(p))
dev.off()




fixdeltansnps=df[df$Delta==1e-6 & df$Nsnps>5000 & df$FPR<0.5,]
theorycopy2=fixdeltansnps
theorycopy2$Power=apply(theorycopy2,1,function(x){relativetheorypower(n=x[1],N=1000,delta=1e-6,alpha=x[2],phi=x[5])})
#write.table(fixdeltansnps,file="emp_relatives.txt",quote=F,row.names=F,col.names=T)
#write.table(theorycopy2,file="theory_relatives.txt",quote=F,row.names=F,col.names=T)
fixdeltansnps$Type<-rep("Empirical",dim(fixdeltansnps)[1])

theorycopy2$Type<-rep("Theory",dim(theorycopy2)[1])
combined=rbind(fixdeltansnps,theorycopy2)



p<-ggplot(fixdeltansnps,aes(x=FPR,y=Power,colour=as.factor(Phi),linetype=Type))+geom_line(lwd=1.5)+facet_wrap(~Nsnps,ncol=3)+theme(strip.text.x = element_text(size=14))+ scale_x_log10(breaks=breakvec,labels=labelvec) +xlab("False Positive Rate")
p<-p+scale_colour_manual(values=cbbPalette,name="Relatedness")
png(paste("E:\\Privacy Code Works\\beacon_distribute\\beacon_distribute\\phi_effect90.png",sep=''),type="cairo",width = 960, height = 960)
print(enhance(p))
dev.off()

p<-ggplot(fixdeltansnps[fixdeltansnps$FPR==0.05,],aes(x=Nsnps,y=Power,linetype=Type,colour=as.factor(Phi)))+geom_line(lwd=1.5)
#p<-p+ scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)) ) 
p<-p+ scale_x_log10(breaks=c(1000,5000,20000,40000,100000),labels=c("1k","5k","20k","40k","100k"))
p<-p+scale_colour_manual(values=cbbPalette,name="Relatedness")+xlab("Number of SNPs")
png(paste("E:\\Privacy Code Works\\beacon_distribute\\beacon_distribute\\phi_effect_power90.png",sep=''),type="cairo",width = 600, height = 480)
print(enhance(p))
dev.off()

#p<-ggplot(combined[combined$Nsnps==5000,],aes(x=FPR,y=Power,linetype=Type,colour=as.factor(Phi)))+geom_line(lwd=1.5) + scale_x_log10(breaks=breakvec,labels=labelvec) +xlab("False Positive Rate")
p<-ggplot(fixdeltansnps[fixdeltansnps$Nsnps==5000,],aes(x=FPR,y=Power,colour=as.factor(Phi)))+geom_line(lwd=1.5) + scale_x_log10(breaks=breakvec,labels=labelvec) +xlab("False Positive Rate")
p<-p+scale_colour_manual(values=cbbPalette,name="Relatedness")
png(paste("E:\\Privacy Code Works\\beacon_distribute\\beacon_distribute\\phi_effect_singlepanel",bias[3],".png",sep=''),type="cairo",width = 600, height = 480)
print(enhance(p))
dev.off()

#ROC curve
fixnsnpsphi=df[df$Phi==1 & df$Nsnps<=40000,]
#  theorycopy3=fixnsnpsphi
#  theorycopy3$Power=apply(theorycopy3,1,function(x){theorypower(n=x[1],N=1000,delta=x[4],alpha=x[2])})
fixnsnpsphi$Type<-rep("Empirical",dim(fixnsnpsphi)[1])
 # theorycopy3$Type<-rep("Theory",dim(theorycopy3)[1])

combined=rbind(fixnsnpsphi)
p<-ggplot(fixnsnpsphi,aes(x=FPR,y=Power,colour=as.factor(Delta)))+geom_line(lwd=1.0)+ylim(0,1)+
  facet_wrap(~Nsnps,ncol=4)+theme(strip.text.x = element_text(size=14))+ scale_x_log10(breaks=breakvec,labels=labelvec) +xlab("False Positive Rate")
p<-p+scale_colour_manual(values=cbbPalette,name="Mismatch Rate")
png(paste("E:\\Privacy Code Works\\beacon_distribute\\beacon_distribute\\delta_effect",bias[6],".png",sep=''),type="cairo",width = 3000, height = 3000,res=200)
print(enhance(p))
dev.off()

p<-ggplot(fixnsnpsphi[fixnsnpsphi$Nsnps==5000,],aes(x=FPR,y=Power,colour=as.factor(Delta)))+geom_line(lwd=1.5)+ scale_x_log10(breaks=breakvec,labels=labelvec) +xlab("False Positive Rate")
p<-p+scale_colour_manual(values=cbbPalette,name="Mismatch Rate")
png(paste("E:\\Privacy Code Works\\beacon_distribute\\beacon_distribute\\delta_effect_singlepanel",bias[6],".png",sep=''),type="cairo",width = 600, height = 480)
print(enhance(p))
dev.off()

