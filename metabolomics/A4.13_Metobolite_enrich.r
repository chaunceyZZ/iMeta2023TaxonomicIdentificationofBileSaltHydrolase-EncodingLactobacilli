arg <- commandArgs(T)
if(length(arg) != 2){
    cat("Argument: pathsum Out_Dir\n")
    quit('no')
}

dataFile = arg[1];
outDir = arg[2];

refFile="/home/Metagenome-2/Ref_database/HMDB/MP_summary.txt";

SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1];
MID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][2];
SID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][3];
VID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][4];
SamID <- paste(SamID,MID,SID,VID,sep=".");

data <- read.table(dataFile,sep="\t")
all.data <- read.table(refFile,sep="\t")

Top <- 12;

pdf(paste(outDir,"/",SamID,".metabo.enriched.pdf",sep=""),width=10,height=7)

Chisq <- matrix(,nrow=nrow(data),ncol=7)
colnames(Chisq) <- c("count1","count2","count3","count4","Pvalue","Adjust_P","-log.P")
rownames(Chisq) <- data[,1]
Chisq[,1] <- data[,2]
Chisq[,2] <- data[,3]
for (i in 1:nrow(data)){
	for (j in 1:nrow(all.data)){
		if (as.character(all.data[j,1])==as.character(data[i,1])){
			temp<-rbind(data[i,2:3],all.data[j,2:3])
			Chisq[i,5] <- signif(chisq.test(temp)$p.value,6)
            Chisq[i,3] <- all.data[j,2]
            Chisq[i,4] <- all.data[j,3]
		}
	}
}
Chisq[,6] <- p.adjust(Chisq[,5],'BH')
Chisq[,7] <- (-log10(abs(Chisq[,6])))
Chisq <- Chisq[Chisq[,6]<0.05,]
write.table(Chisq,paste(outDir,"/",SamID,".Chisq.txt",sep=""),quote = F,sep ="\t",row.names=T,col.names=F)
Chisq <- Chisq[order(Chisq[,7],decreasing=T),]
Chisq <- cbind(Chisq[1:Top,1],Chisq[1:Top,5:7])
Chisq <- Chisq[order(Chisq[,1],decreasing=F),]
cPal <- colorRampPalette(c("blue", "white","red" ))
mycolors <- cPal(1.2*nrow(Chisq))[as.numeric(cut(Chisq[,4],breaks = 1.2*nrow(Chisq)))]
par(mar=c(5,20,5,3))
bar.pos<-barplot(Chisq[,1],space=0.5,horiz=T,col=mycolors,las=2,axes=F,cex.names=0.9,xlim=c(0,(max(Chisq[,1])*1.1)),border=NA,axisnames=T)
axis(2,at=apply(matrix(c(0,as.vector(bar.pos),max(bar.pos)+3),nrow=1),2,median),labels=NA,tck=-0.02)
axis(3,at=c(-4,seq(0,(max(Chisq[,1]*1.1)),0.05)),outer=F,line=-1,xpd=T,tick=T,labels=F,tck=0)
axis(3,at=seq(0,(max(Chisq[,1]*1.1)),2),outer=F,line=-1,xpd=T,tck=-0.008,lwd=0,lwd.tick=1,cex.axis=0.8)
mtext('Enriched pathways',side=3,line=-0.8,adj=1,cex=1.2,at=-1)
#mtext("Total count metabolites",side=3,cex=1,line=-1.2)

par(new=T)
par(new=T,fig=c(0.7,0.9,0.05,0.15),mar=c(0.1,1,2,1))
image(as.matrix(seq(min(Chisq[,4]),max(Chisq[,4]),length.out=500)),axes=F,col=cPal(500))       #图例
mtext(expression('-Log'[10]*"(adjusted_P)"),side=3,line=0,cex=1)
max <- round(max(Chisq[,4]))
min <- round(min(Chisq[,4]))
mtext(seq(min,max,1),at=seq(0,1,(1/(max-min))),side=1,line=0,cex=1)
dev.off()


