arg <- commandArgs(T)
if(length(arg) != 2){
        cat("Argument: Data_File Out_Dir\n")
        quit('no')
}

dataFile <- arg[1];
outDir <- arg[2];

SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1];
library(plotrix)

OTU<-read.table(dataFile,header=T);
OTU.count<-apply(OTU[,-1]>0,1,sum);
OTU.frq<-OTU.count/ncol(OTU[,-1]);
hist(OTU.frq,plot=F);
OTU.frq<-OTU.count/ncol(OTU[,-1])
OTU.reads<-apply(OTU[,-1],1,sum)
cutoff<-seq(0,1,0.05)
OTU.coreDist<-matrix(nrow=length(cutoff),ncol=2)
for(i in 1:length(cutoff)){
	OTU.coreDist[i,1] = sum(OTU.frq>=cutoff[i])
	OTU.coreDist[i,2] = sum(OTU.reads[OTU.frq>=cutoff[i]])
}
OTU.coreDist.frq=OTU.coreDist
OTU.coreDist.frq[,1]=OTU.coreDist[,1]/OTU.coreDist[1,1]
OTU.coreDist.frq[,2]=OTU.coreDist[,2]/OTU.coreDist[1,2]


y.tick=0.002
OTU.coreDist.4plot.otu=OTU.coreDist.frq
OTU.coreDist.4plot.otu[,2]=0
OTU.coreDist.4plot.reads=OTU.coreDist.frq
OTU.coreDist.4plot.reads[,1]=0
OTU.coreDist.4plot.otu[,1][OTU.coreDist.frq[,1]>0.02]<-(0.01+y.tick*1)+(OTU.coreDist.4plot.otu[,1][OTU.coreDist.frq[,1]>0.02]-(0.01+y.tick*1))/(1-(0.01+y.tick*1))*(y.tick*1)
OTU.coreDist.4plot.otu[,1][OTU.coreDist.frq[,1]<=0.02 & OTU.coreDist.frq[,1]>0.01]<-0.01+(OTU.coreDist.4plot.otu[,1][OTU.coreDist.frq[,1]<=0.02 & OTU.coreDist.frq[,1]>0.01]-0.01)/(0.02-0.01)*(y.tick*1)

pdf(paste(outDir,"/",SamID,".CommonOTUs.pdf",sep=""),width=8,height=5)
barplot(t(OTU.coreDist.4plot.otu),beside=T,axes=F,col=c('darkseagreen4','darkslategray'))
axis(1,at=barplot(t(OTU.coreDist.4plot.otu),beside=T,axes=F,plot=F)[1,],labels=paste(">=",seq(0,100,5),"%",sep=''),cex.axis=0.65)
axis(2,at=seq(0,0.014,0.002),labels=c(seq(0,0.01,0.002),0.02,1),tck=-0.01,col='darkseagreen4')
axis.break(2,breakpos=0.01,style='zigzag',breakcol='darkseagreen4')
axis.break(2,breakpos=0.012,style='zigzag',breakcol='darkseagreen4')
par(new=T)
barplot(t(OTU.coreDist.4plot.reads),beside=T,axes=F,col=c('darkseagreen4','darkslategray'))
legend(x='topright',bty='n',fill=c('darkseagreen4','darkslategray'),c('Common OTUs','Reads in OTUs'),cex=0.85)
axis(4,at=seq(0,1,0.1),labels=seq(0,1,0.1),tck=-0.01,col='darkslategray')
dev.off()
