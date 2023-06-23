arg <- commandArgs(T)
if(length(arg) != 5){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}

print(arg)
#arg <-c( "/home/wmj/02_Metagenome/Metabolome/20190904_SZW/Temp/1.2_raw_Distance/FHC.spearman.txt", "/home/wmj/02_Metagenome/Metabolome/20190904_SZW/Temp/bin/FHC-grouping.info" , "/home/wmj/02_Metagenome/Metabolome/20190904_SZW/Temp/1.4_raw_MDS"  , "2", "/home/wmj/02_Metagenome/Metabolome/20190904_SZW/Temp/bin/FHC-color.txt" )
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

grpInfo <- read.table(groupFile,header=T);

FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- FileName[1];
MID <- FileName[2];

Beta <- as.matrix(data.frame(read.table(dataFile,header=T,row.names = 1)));

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

Groups <- as.character(levels(as.factor(grpInfo[,grp_col])))
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


mds <- isoMDS(Beta,k=2);

pdf(paste(outDir,"/",SamID,".",MID,".mds.pdf",sep=""),width=4,height=4);

plot(mds$points[,1],mds$points[,2],col = Color.Group,pch=Symbol.Group,xlab="Dimension1",ylab="Dimension2", main=paste("Stress=",sprintf("%.3f",mds$stress),sep=""),cex=0.9);

s.class(mds$points, fac=as.factor(groupname),grid=F, addaxes=F,axesell =T,label=Groups,col=unique(Color.Group),pch=Symbol.Group,add.plot=T);

dev.off()
