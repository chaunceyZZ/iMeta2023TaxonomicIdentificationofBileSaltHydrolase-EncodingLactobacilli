arg <- commandArgs(T)
if(length(arg) < 4){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}
print(arg)
#arg <-c("/home/wmj/02_Metagenome/16s/HTQ_190723/14_Statistics/beta/Total/spearman_approx_Total.gut.sub2W2.filter_4_noChimera_rmAlignFail.otu_table_even15000.txt", "/home/wmj/02_Metagenome/16s/HTQ_190723/bin/Total-color.txt" , "/home/wmj/02_Metagenome/16s/HTQ_190723/Analysis/8_OTUs_beta/Total"  ,"2" )

library(vegan);
library(MASS);
#library(cluster);
#library(clusterSim);
library(ade4);
#library(fpc);


dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);

#dataFile <- "abund_jaccard.gut.resample.otu.txt";
#groupFile <- "gut-grouping.info";
#outDir <- "./";
#grp_col <- 2;


grpInfo <- read.table(groupFile,header=T,check.names = FALSE, comment.char = "")
rownames(grpInfo) <- as.character(grpInfo[,1])
FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
MetID <- FileName[1];


Beta <- as.matrix(data.frame(read.table(dataFile,header=T,row.names = 1,check.names = FALSE)));
ii <- intersect(rownames(grpInfo),rownames(Beta))
if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),ii)
	print(ii1)
}


Beta <- as.matrix(Beta[ii,ii])
grpInfo <- grpInfo[ii,]

Beta <- Beta[order(rownames(Beta),decreasing=F),];
Beta <- Beta[,order(colnames(Beta),decreasing=F)];

groupname <- c();
for(i in 1:length(rownames(Beta)))
{
	groupname <- append(groupname,as.character(grpInfo[grep(paste("^",rownames(Beta)[i],"$",sep=""),grpInfo[,1]),grp_col]));
}

rownames(Beta) <- groupname;
colnames(Beta) <- groupname;

Groups <- as.character(unique(grpInfo[,grp_col]))
G_cc=as.character(unique(grpInfo[,3]))
names(G_cc)=Groups
#gg <- read.table(arg[5],header=T,check.names = FALSE)
#rownames(gg) <- as.character(gg[,1])
#ii <- intersect(rownames(grpInfo),rownames(gg))
#gg <- gg[ii,]
gg <- grpInfo

Groups.1 <-as.character(unique(gg[,2])) 


Beta.Group <- c();
Symbol.Group <- c();
Color.Group <- c();

for(i in 1:length(Groups))
{
	Beta.Group[[i]] <- Beta[grep(paste("^",Groups[i],"$",sep=""),rownames(Beta)),grep(paste("$",Groups[i],"$",sep=""),colnames(Beta))];
	Symbol.Group <- append(Symbol.Group,as.vector(rep((17+i),nrow(Beta.Group[[i]]))));
	Color.Group <- append(Color.Group,as.vector(rep(rainbow(length(Groups))[i],nrow(Beta.Group[[i]]))));
}

mds <- isoMDS(Beta,k=2);

pdf(paste(outDir,"/",MetID,".resample.otu.mds.pdf",sep=""),width=4,height=4);
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=G_cc
plot(mds$points[,1],mds$points[,2],col = col.line[Groups],pch=Symbol.Group,xlab="Dimension1",ylab="Dimension2", main=paste("Stress=",sprintf("%.3f",mds$stress),sep=""),cex=0.9);

s.class(mds$points, fac=as.factor(groupname),grid=F, addaxes=F,axesell =T,label=Groups.1,col=col.line[Groups],pch=Symbol.Group,add.plot=T);

dev.off()
