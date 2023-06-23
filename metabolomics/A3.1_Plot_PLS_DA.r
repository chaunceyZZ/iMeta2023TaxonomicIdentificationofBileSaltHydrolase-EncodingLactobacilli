arg <- commandArgs(T)
if(length(arg) != 5){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}
print(arg)
#arg <-c("/home/wmj/Metagenome-2/Metabolome/HFQ/HFQ/HFQ_New-1/P_Dis/0.0_raw_data/FBM.data.txt", "/home/wmj/Metagenome-2/Metabolome/HFQ/HFQ/HFQ_New-1/P_Dis/bin/FBM-grouping.info" , "/home/wmj/Metagenome-2/Metabolome/HFQ/HFQ/HFQ_New-1/P_Dis/3.1_PLS_DA" , "2", "/home/wmj/Metagenome-2/Metabolome/HFQ/HFQ/HFQ_New-1/P_Dis/bin/Fall-color.txt"    )

library(mixOmics);
library(ade4);
library(vegan);
library("car")

dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);

#dataFile <- "total.data.txt";
#groupFile <- "total-grouping.info";
#outDir <- "./";
#grp_col <- 2;


grpInfo <- read.table(groupFile,header=T);
FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- FileName[1];
Data <- read.table(dataFile,header=T,row.names = 1);
Data.tm <- t(as.matrix(Data)); 
groupname <- c();

for(i in 1:length(colnames(Data)))
{
	groupname <- append(groupname,as.character(grpInfo[grep(paste("^",colnames(Data)[i],"$",sep=""),as.character(grpInfo[,1])),grp_col]));
}

Data.plsda <- plsda(Data.tm ,groupname,ncomp=2)

pdf(paste(outDir,"/",SamID,".pls-da.pdf",sep=""),height=8,width=8);

plotIndiv(Data.plsda,ind.names = F, ellipse = T, centroid= F, title= "PLS-DA Plot", legend =T,style="ggplot2")

dev.off()


pdf(paste(outDir,"/",SamID,".point.pls-da.pdf",sep=""),height=8,width=8);
xy <- Data.plsda$variates
xx <- xy[[1]][,1]
yy <- xy[[1]][,2]

lab <- Data.plsda$explained_variance
x.lab <-round(lab[[1]][1]*100,0)
y.lab <-round(lab[[1]][2]*100,0)
col1 <- read.table(arg[5],header=T,row.names=1,sep="\t",comment.char = "")
ii <- names(xx)
col <-as.character(col1[ii,2])



par(mgp=c(2.5,1,0),mar=c(10,5,2,2))
plot(x=xx,y=yy,pch=rep(20,length(col)),col=col,main=paste(SamID," PLS-DA Plot",sep=""),xlab=paste("X−variate 1:",x.lab,"% expl. var",sep=""),ylab=paste("X−variate 2:",y.lab,"% expl. var",sep=""),cex.lab=2,cex.main=2,cex.axis=2,xlim=c(-50,35),ylim=c(-30,30))


dataEllipse(x=as.numeric(xx),y=as.numeric(yy),groups=factor(groupname,levels=unique(groupname),labels=unique(groupname)),levels = c(0.95), add=TRUE, col =unique(col) , lwd = 1,plot.points=FALSE,fill=TRUE,center.cex=0.2,group.labels="",fill.alpha = 0)

legend("topleft",pch=rep(20,length(unique(col))),col=unique(col),bty="n",legend=unique(groupname),cex=1.8,pt.cex=1.8)


dev.off()






