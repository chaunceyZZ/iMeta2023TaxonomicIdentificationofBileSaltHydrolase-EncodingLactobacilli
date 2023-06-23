arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}
print(arg)
#arg <-c("/home/Project/16S/LSZ_20190908/Analysis/12_Sig_Distance/Total/Total.Utest.oneSig.spearman.txt", "/home/Project/16S/LSZ_20190908/bin/Total-color.txt" , "/home/Project/16S/LSZ_20190908/Analysis/14_Sig_MDS/Total"   , "2"  )
dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);

#dataFile <- "total.Utest.spearman.txt";
#groupFile <- "total-grouping.info";
#outDir <- "./";
#grp_col <- 2;

library(vegan);
library(MASS);
#library(cluster);
#library(clusterSim);
library(ade4);
#library(fpc);

grpInfo <- read.table(groupFile,header=T,check.names = FALSE, comment.char = "")
rownames(grpInfo)=as.character(grpInfo[,1])
FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- paste(FileName[1],FileName[2],sep=".");
MID <- FileName[3];
LID <- FileName[4];
DID <- FileName[5];
GID <- colnames(grpInfo)[grp_col];

Beta <- as.matrix(data.frame(read.table(dataFile,header=T,row.names = 1,check.names = FALSE),check.names = FALSE),check.names = FALSE);


ii <- intersect(rownames(grpInfo),colnames(Beta))
if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),ii)
	print(ii1)
}
grpInfo=grpInfo[ii,]
Beta=Beta[ii,ii]



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

col.line=unique(as.character(grpInfo[,3]))
names(col.line)=unique(as.character(grpInfo[,2]))

for(i in 1:length(Groups))
{
	Beta.Group[[i]] <- Beta[grep(paste("^",Groups[i],"$",sep=""),rownames(Beta)),grep(paste("^",Groups[i],"$",sep=""),colnames(Beta))];
	Symbol.Group <- append(Symbol.Group,as.vector(rep((17+i),nrow(Beta.Group[[i]]))));
#	Color.Group <- append(Color.Group,as.vector(rep(rainbow(length(Groups))[i],nrow(Beta.Group[[i]]))));
	Color.Group <- append(Color.Group,as.vector(rep(col.line[i],nrow(Beta.Group[[i]]))));
}

mds <- isoMDS(Beta,k=2);



pdf(paste(outDir,"/",SamID,".",MID,".",LID,".",DID,".mds.pdf",sep=""),width=4,height=4);

plot(mds$points[,1],mds$points[,2],col = Color.Group,pch=Symbol.Group,xlab="Dimension1",ylab="Dimension2", main=paste("Stress=",sprintf("%.3f",mds$stress),sep=""),cex=0.9);

#s.class(mds$points, fac=as.factor(groupname),grid=F, addaxes=F,axesell =T,label=Groups,col=rainbow(length(Groups)),pch=Symbol.Group,add.plot=T);
s.class(mds$points, fac=as.factor(groupname),grid=F, addaxes=F,axesell =T,label=Groups,col=col.line[Groups],pch=Symbol.Group,add.plot=T);

dev.off()
