library(ggplot2);
arg <- commandArgs(T)
if(length(arg) != 4){
	cat("Argument: dataFile GroupInfo Out_Dir Group_column Y C\n")
	quit('no')
}
print(arg)
dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);
#fc <- arg[5:6];

#dataFile <- "high.data.txt";
#groupFile <- "high-grouping.info";
#outDir <- "./";
#grp_col <- 2

data <- read.table(dataFile,header=T,row.names=1)
grp <- read.table(groupFile,header=T);
fc <- as.character(unique(grp[,grp_col]))
SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1];
MID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][2];
SID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][3];
VID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][4];
SamID <- paste(SamID,MID,SID,VID,sep=".");

stat <- matrix(,nrow=nrow(data),ncol=4)
colnames(stat) <- c("fold change","pvalue","logFC","-logP")
rownames(stat) <- rownames(data)

groupname <- c();
for(i in 1:length(colnames(data)))
{
	groupname <- append(groupname,as.character(grp[grep(colnames(data)[i],as.character(grp[,1])),grp_col]));
}

colnames(data) <- groupname
Groups <- as.character(levels(as.factor(grp[,grp_col])))

if(length(Groups) >2)
{
	quit('no');
}

data.Group <- list()
for(i in 1:length(Groups))
{
	data.Group[[i]] <- data[,grep(fc[i],colnames(data))];
}


for (i in 1:nrow(data)){
	stat[i,1] <- mean(as.numeric(data.Group[[1]][i,]))/mean(as.numeric(data.Group[[2]][i,]))
	stat[i,2] <- wilcox.test(as.numeric(data.Group[[1]][i,]),as.numeric(data.Group[[2]][i,]))$p.value
}
stat[,3] <- (log2(abs(stat[,1])))
stat[,4] <- (-log10(abs(stat[,2])))
write.csv(stat,paste(outDir,"/",SamID,".Volcano.csv",sep=""))
stat<-data.frame(stat);
pdf(paste(outDir,"/",SamID,".Volcano.pdf",sep=""),w=8,h=8);
threshold <- as.factor((stat$logFC>1 | stat$logFC<(-1)) & stat$pvalue<0.05)
par(mar=c(4,4,2,2),mgp=c(2,1,0))
mycol <- rep('grey',length=nrow(stat))
mycol[grep('TRUE',threshold)] <- 'red'
plot(x=stat[,3],y=stat[,4],col=mycol,main=paste0(fc[1]," vs ",fc[2]),xlab="log2(Fold Change)",ylab="-log10(p-value)",bty='l',pch=20)
dev.off()
