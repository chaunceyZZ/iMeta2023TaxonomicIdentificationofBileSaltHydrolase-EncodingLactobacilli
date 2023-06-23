arg <- commandArgs(T)
if(length(arg) != 5){
        cat("Argument: Data_File Group_File Out_Dir Group_column Color\n")
        quit('no')
}
print(arg)
#arg <-c("/home/Project/Metabolome/HFQ/HFQ_New-1/P_Dis/4.3_VIP_Sig_Dis/FDN.Utest.Sig.pearson.VIP.txt","/home/Project/Metabolome/HFQ/HFQ_New-1/P_Dis/bin/FDN-grouping.info", "/home/Project/Metabolome/HFQ/HFQ_New-1/P_Dis/4.4_VIP_Sig_PCA" , "2"  ,"/home/Project/Metabolome/HFQ/HFQ_New-1/P_Dis/bin/FDN-color.txt" )

dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);

#dataFile <- "abund_jaccard_Total.resample.otu.txt";
#groupFile <- "gut-grouping.info";
#outDir <- "./";
#grp_col <- 2;

#library(vegan);
#library(MASS);
#library(cluster);
#library(clusterSim);
library(ade4);
library(fpc);
library(ggplot2);
library("car")

grpInfo <- read.table(groupFile,header=T);

FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- FileName[1];
MetID <- FileName[2];
LID <- FileName[3];
DID <- FileName[4];

Beta <- read.table(dataFile,header=T,sep="\t")
nn <- as.character(Beta[,1])
dd <- as.matrix(Beta[,-1])
rownames(dd) <- nn

Beta <- dd

Beta <- Beta[order(rownames(Beta),decreasing=F),];
Beta <- Beta[,order(colnames(Beta),decreasing=F)];

groupname <- c();
ii <-c()
for(i in 1:length(rownames(Beta)))
{
	groupname <- append(groupname,as.character(grpInfo[grep(paste("^",rownames(Beta)[i],"$",sep=""),as.character(grpInfo[,1])),grp_col]));
	ii <-c(ii,as.character(grpInfo[grep(paste("^",rownames(Beta)[i],"$",sep=""),as.character(grpInfo[,1])),1]))
}

rownames(Beta) <- groupname;
colnames(Beta) <- groupname;

Groups <- unique(groupname)
Beta.Group <- c();
Symbol.Group <- c();
Color.Group <- c();

p.info <- read.table(arg[5],header=T,row.names=1,comment.char = "")
p.info <- p.info[ii,]
cc=unique(p.info)
c_col=as.character(cc[,2])
names(c_col)=as.character(cc[,1])
for(i in 1:length(Groups))
{
	Beta.Group[[i]] <- Beta[grep(paste("^",Groups[i],"$",sep=""),rownames(Beta)),grep(paste("^",Groups[i],"$",sep=""),colnames(Beta))];
	Symbol.Group <- append(Symbol.Group,as.vector(rep((17+i),nrow(Beta.Group[[i]]))));
	Color.Group <- append(Color.Group,rep(c_col[Groups[i]],nrow(Beta.Group[[i]])));
}


pc <- princomp(Beta);
#write.csv(pc$loadings[,1:2],paste(outDir,"/",SamID,".",MetID,".",LID,".",DID,".VIP.pca.csv",sep=""),row.names=T)
pdf(paste(outDir,"/",SamID,".",MetID,".",LID,".",DID,".VIP.pca.pdf",sep=""),height=8,width=8);

plot(pc$loadings[,1],pc$loadings[,2],col =Color.Group,pch=Symbol.Group,xlab=paste("PC1(",sprintf("%.1f",pc$sde[1]^2/sum(pc$sde^2)*100),"%)",spe=""),ylab=paste("PC2(",sprintf("%.1f",pc$sde[2]^2/sum(pc$sde^2)*100),"%)",spe=""),cex=2);

s.class(pc$loadings[,1:2], fac=factor(groupname,levels=unique(groupname),labels=unique(groupname)),grid=F, addaxes=F,axesell =F,label=Groups,col=unique(Color.Group),pch=Symbol.Group,add.plot=T);

dev.off()


### PCA.AERA
pdf(paste(outDir,"/",SamID,".",MetID,".",LID,".",DID,".VIP.pca.area.pdf",sep=""),height=8,width=8);
par(mgp=c(2.5,1,0),mar=c(10,5,2,2),tcl=-0.2)
test <- data.frame(pc$loading[,1:2],group=groupname,col=ii,lab=ii)
p.info <- read.table(arg[5],header=T,row.names=1,comment.char = "")
p.info <- p.info[ii,]
#p.info[,2] <- paste("#",p.info[,2],sep="")

test[,4] <- p.info[,2]
test[,5] <- p.info[,1]

x <- test[,1]
y <- test[,2]
plot(x=x,y=y,pch=rep(20,nrow(test)),col=as.character(test[,4]),xlab=paste("PC1(",sprintf("%.1f",pc$sde[1]^2/sum(pc$sde^2)*100),"%)",spe=""),ylab=paste("PC2(",sprintf("%.1f",pc$sde[2]^2/sum(pc$sde^2)*100),"%)",spe=""),cex.lab=2,cex.main=2,cex.axis=2,xlim=c(-0.2,0.2),ylim=c(-0.3,0.3),cex=2)



dataEllipse(x=as.numeric(x),y=as.numeric(y),groups=factor(test[,3],levels=unique(test[,3]),labels=unique(test[,3])),levels = c(0.95), add=TRUE, col =as.character(unique(test[,4])) , lwd = 1,plot.points=FALSE,fill=TRUE,center.cex=0.2,group.labels="",fill.alpha = 0.3)
#abline(h=0,lwd=2)
#abline(v=0,lwd=2)
legend("bottom",pch=rep(20,length(unique(test[,4]))),col=as.character(unique(test[,4])),bty="n",legend=as.character(unlist(unique(test[,5]))),cex=1.8,pt.cex=1.8,xpd=T,inset=-0.3,horiz=T)

dev.off()


### PCA.point
pdf(paste(outDir,"/",SamID,".",MetID,".",LID,".",DID,".VIP.pca.point.pdf",sep=""),height=8,width=8);
par(mgp=c(2.5,1,0),mar=c(10,5,2,2),tcl=-0.2)
test <- data.frame(pc$loading[,1:2],group=groupname,col=ii,lab=ii)
p.info <- read.table(arg[5],header=T,row.names=1,comment.char = "")
p.info <- p.info[ii,]
#p.info[,2] <- paste("#",p.info[,2],sep="")

test[,4] <- p.info[,2]
test[,5] <- p.info[,1]

x <- test[,1]
y <- test[,2]
plot(x=x,y=y,pch=rep(20,nrow(test)),col=as.character(test[,4]),xlab=paste("PC1(",sprintf("%.1f",pc$sde[1]^2/sum(pc$sde^2)*100),"%)",spe=""),ylab=paste("PC2(",sprintf("%.1f",pc$sde[2]^2/sum(pc$sde^2)*100),"%)",spe=""),cex.lab=2,cex.main=2,cex.axis=2,xlim=c(-0.3,0.3),ylim=c(-0.4,0.4),cex=2)



dataEllipse(x=as.numeric(x),y=as.numeric(y),groups=factor(test[,3],levels=unique(test[,3]),labels=unique(test[,3])),levels = c(0.95), add=TRUE, col =as.character(unique(test[,4])) , lwd = 1,plot.points=FALSE,fill=TRUE,center.cex=0.2,group.labels="",fill.alpha = 0)
#abline(h=0,lwd=2)
#abline(v=0,lwd=2)
legend("bottom",pch=rep(20,length(unique(test[,4]))),col=as.character(unique(test[,4])),bty="n",legend=as.character(unlist(unique(test[,5]))),cex=1.8,pt.cex=1.8,xpd=T,inset=-0.3,horiz=T)

dev.off()


