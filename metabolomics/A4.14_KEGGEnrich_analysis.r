
arg <- commandArgs(T)
if(length(arg) != 4){
	cat("Argument: dataFile GroupInfo Out_Dir Group_column Y C\n")
	quit('no')
}

print(arg)
#arg <-c("/project/Metabolome/HFQ/HFQ_New-1/P_Dis/ID_Metabolites/FDN.ID.txt","/project/Metabolome/HFQ/HFQ_New-1/P_Dis/4.1_VIP_Sig/FDN.Utest.Sig.VIP.txt","/project/Metabolome/HFQ/HFQ_New-1/P_Dis/bin/FDN-grouping.info" ,"/project/Metabolome/HFQ/HFQ_New-1/P_Dis/4.10_ID_KEGG")
library("clusterProfiler")
library("DOSE")
library("org.Hs.eg.db")   ### "org.Mm.eg.db" for Mouse, if for Human should be "org.Hs.eg.db"
library("pathview")
library(ggplot2);
###


setwd(arg[4])
id <- read.table(arg[1],header=T,sep='\t')
dir3 <- dirname(arg[3])
ref <- read.table(paste0(dir3,"/",'ALL_HMDB_gene.txt'),sep='\t')
n.name <- strsplit(basename(arg[1]),".",fixed=T)[[1]][1]
gene <- ref[ref[,1] %in% id[,2],] #ref是否包含在id中
#id[id[,3] %in% ref[,1],]

abun <- read.table(arg[2],header=T,row.names=1)
abun <- abun[(as.character(id[id[,2] %in% ref[,1],1])),]
flag <- c()
gg <- read.table(arg[3],header=T,row.names=1)
Groups <- as.character(unique(gg[,1]))
for (i in 1:nrow(abun))
{
	i1 <- rownames(gg)[which(as.character(gg[,1])==Groups[1])]
	i2 <- rownames(gg)[which(as.character(gg[,1])==Groups[2])]
	t1 <- mean(as.numeric(abun[i,i1]))
	t2 <- mean(as.numeric(abun[i,i2]))
  	if(t1 >= t2)
   	{
       		 flag[i] <- Groups[1]
    	}else{
       	 	flag[i] <- Groups[2]
    	}
}
dat <- cbind(id[id[,2] %in% ref[,1],],flag=flag)

flag2 <- c()
for (i in 1:nrow(gene))
{
	flag2[i] <- as.character(dat[grep(gene[i,1],dat[,2]),4][1])
}

data <- data.frame(cbind(flag2,as.character(gene[,3])))
data <- data[unique(data[,2]),]
gene <- as.character(data[,2]);
gene.df <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Hs.eg.db)
gene.id <- gene.df[,2];

gene.gr <- c();
for(i in 1:length(gene.df[,1]))
{
	gene.gr <- append(gene.gr,as.character(data[grep(gene.df[i,1],data[,2])[1],1]));
}

NewData <- cbind(gene.gr,gene.id);
#Groups1 <- as.character(levels(as.factor(NewData[,1])));
Groups1 <- as.character(unique(NewData[,1]));
for(i in 1:length(Groups1))
{
	gene.sub <- NewData[grep(as.character(Groups1[i]),as.character(NewData[,1])),2];	
#	gene.sub <- NewData[grep(Groups1[i],NewData[,1]),2];
	ekk <- enrichKEGG(gene=gene.sub,organism="hsa",pvalueCutoff=0.05,pAdjustMethod = "BH");
	if(is.null(ekk)) next
	else{
	write.csv(summary(ekk),paste(n.name,".",Groups1[i],".VIP.KEGGenrich.csv",sep=""),row.names=F);
	pdf(paste(n.name,".",Groups1[i],".VIP.KEGGcnetplot.pdf",sep=""),h=15,w=15);
       	cnetplot(ekk,categorySize="pvalue");
	dev.off();
	}
}

pdf(paste(n.name,".KEGGenriched.barplot.pdf",sep=""),width=8,height=5);
ekk <- enrichKEGG(gene=gene.id,organism="hsa",pvalueCutoff=0.05,pAdjustMethod = "BH");
write.csv(summary(ekk),paste(n.name,".VIP.barplot.KEGGenrich.csv",sep=""),row.names=F);
barplot(ekk,drop=TRUE,showCategory = 12);
dev.off()

mydf <- data.frame(ID = gene.id, Group = gene.gr);

gene.comp <- compareCluster(ID~Group,data=mydf,fun="enrichKEGG",organism="hsa",pvalueCutoff=0.05,pAdjustMethod = "BH")
pdf(paste(n.name,".KEGGenriched.dot.pdf",sep=""),width=8,height=5)
plot(gene.comp);
dev.off()




