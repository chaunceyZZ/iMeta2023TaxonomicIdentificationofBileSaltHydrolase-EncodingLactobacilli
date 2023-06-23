library("clusterProfiler")
library("DOSE")
library("org.Hs.eg.db")   ### "org.Mm.eg.db" for Mouse, if for Human should be "org.Hs.eg.db"
library("pathview")
arg <- commandArgs(T)
if(length(arg) != 2){
        cat("Argument: Data_File Out_Dir\n")
        quit('no')
}

dataFile <- arg[1];
outDir <- arg[2];

#dataFile <- "FM.Utest.Sig.VIP.enriched.gene.nred.2.dat";
#outDir  <- "./";

FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- FileName[1];
MID <- FileName[2];
LID <- FileName[3];


data <- read.table(dataFile,header=F);
gene <- as.character(data[,2]);
gene.df <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Hs.eg.db)
gene.id <- gene.df[,2];

gene.gr <- c();
for(i in 1:length(gene.df[,1]))
{
	for(j in 1:length(data[,2]))
	{
		if(gene.df[i,1] == data[j,2])
		{
			gene.gr <- append(gene.gr,as.character(data[j,1]));
		}
	}
}

NewData <- cbind(gene.gr,gene.id);
Groups <- as.character(levels(as.factor(NewData[,1])));
for(i in 1:length(Groups))
{
	gene.sub <- NewData[grep(Groups[i],NewData[,1]),2];
	ekk <- enrichKEGG(gene=gene.sub,organism="hsa",pvalueCutoff=0.05,pAdjustMethod = "BH");
	if(is.null(ekk)) next
	else{
	write.csv(summary(ekk),paste(outDir,"/",SamID,".",MID,".",LID,".",Groups[i],".VIP.KEGGenrich.csv",sep=""),row.names=F);
	pdf(paste(outDir,"/",SamID,".",MID,".",LID,".",Groups[i],".VIP.cnetplot.pdf",sep=""),h=15,w=15);
        cnetplot(ekk,categorySize="pvalue");
        dev.off();
	for(j in 1:1)
        {
                tt <- pathview(gene.data=gene.sub,pathway.id=ekk$ID[j],species="hsa", out.suffix=paste(SamID,".",MID,".",LID,".",Groups[i],sep=""));
        }
	}
}

pdf(paste(outDir,"/",SamID,".enriched.barplot.pdf",sep=""),width=8,height=5);
ekk <- enrichKEGG(gene=gene.id,organism="hsa",pvalueCutoff=0.05,pAdjustMethod = "BH");
barplot(ekk,drop=TRUE,showCategory = 12);

mydf <- data.frame(ID = gene.id, Group = gene.gr);

gene.comp <- compareCluster(ID~Group,data=mydf,fun="enrichKEGG",organism="hsa",pvalueCutoff=0.05,pAdjustMethod = "BH")
pdf(paste(outDir,"/",SamID,".",MID,".",LID,".VIP.enriched.dot.pdf",sep=""),width=8,height=5)
plot(gene.comp);


dev.off()
