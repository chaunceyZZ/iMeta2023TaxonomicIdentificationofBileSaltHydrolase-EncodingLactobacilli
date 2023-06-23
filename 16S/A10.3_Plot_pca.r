arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}
print(arg)
#arg <-c("/home/wmj/Metagenome-2/16s/HK/AB/Analysis/12_Sig_Distance/Group/All.allSig.data.bray-curtis.txt", "/home/wmj/Metagenome-2/16s/HK/AB/Analysis/4_common_OTUs/Group/All.group.txt","/home/wmj/Metagenome-2/16s/HK/AB/Analysis/13_Sig_PCA/Group", "2" )

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

grpInfo <- read.table(groupFile,header=T,check.names = FALSE, comment.char = "")
rownames(grpInfo)=as.character(grpInfo[,1])
FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- paste(FileName[1],FileName[2],sep=".");
MetID <- FileName[3];
LID <- FileName[4];
DID <- FileName[5];

Beta <- as.matrix(data.frame(read.table(dataFile,header=T,row.names = 1,check.names = FALSE),check.names = FALSE),check.names = FALSE);
ii <- intersect(as.character(grpInfo[,1]),rownames(Beta))
if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),ii)
	print(ii1)
}
grpInfo=grpInfo[ii,]
Beta<-Beta[ii,ii]

Beta <- Beta[order(rownames(Beta),decreasing=F),];
Beta <- Beta[,order(colnames(Beta),decreasing=F)];

groupname <- c();
for(i in 1:length(rownames(Beta)))
{
	groupname <- append(groupname,as.character(grpInfo[grep(paste("^",rownames(Beta)[i],"$",sep=""),grpInfo[,1]),grp_col]));
}

rownames(Beta) <- groupname;
colnames(Beta) <- groupname;

#Groups <- as.character(levels(as.factor(grpInfo[,grp_col])))
Groups <- as.character(unique(grpInfo[,grp_col]))
Beta.Group <- c();
Symbol.Group <- c();
Color.Group <- c();
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=unique(as.character(grpInfo[,3]))
names(col.line)=unique(as.character(grpInfo[,2]))

for(i in 1:length(Groups))
{
	Beta.Group[[i]] <- Beta[grep(paste("^",Groups[i],"$",sep=""),rownames(Beta)),grep(paste("^",Groups[i],"$",sep=""),colnames(Beta))];
	Symbol.Group <- append(Symbol.Group,as.vector(rep((17+i),nrow(Beta.Group[[i]]))));
	Color.Group <- append(Color.Group,as.vector(rep(col.line[i],nrow(Beta.Group[[i]]))));
}

pc <- princomp(Beta);

pdf(paste(outDir,"/",SamID,".",MetID,".",LID,".",DID,".pca.pdf",sep=""),width=4,height=4);

plot(pc$loadings[,1],pc$loadings[,2],col = Color.Group,pch=Symbol.Group,xlab=paste("PC1(",sprintf("%.1f",pc$sde[1]^2/sum(pc$sde^2)*100),"%)",spe=""),ylab=paste("PC2(",sprintf("%.1f",pc$sde[2]^2/sum(pc$sde^2)*100),"%)",spe=""),cex=0.9);

s.class(pc$loadings[,1:2], fac=as.factor(groupname),grid=F, addaxes=F,axesell =T,label=Groups,col=col.line[Groups],pch=Symbol.Group,add.plot=T);

dev.off()

pdf(paste(outDir,"/",SamID,".",MetID,".",LID,".",DID,".pca.nocicle.pdf",sep=""),width=4,height=4);

plot(pc$loadings[,1],pc$loadings[,2],col = Color.Group,pch=Symbol.Group,xlab=paste("PC1(",sprintf("%.1f",pc$sde[1]^2/sum(pc$sde^2)*100),"%)",spe=""),ylab=paste("PC2(",sprintf("%.1f",pc$sde[2]^2/sum(pc$sde^2)*100),"%)",spe=""),cex=0.9);


dev.off()

