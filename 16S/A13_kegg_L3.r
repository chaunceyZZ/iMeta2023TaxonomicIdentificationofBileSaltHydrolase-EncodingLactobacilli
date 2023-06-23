arg <- commandArgs(T)
if(length(arg) < 1){
	cat("Argument: Out_Dir \n")
	quit('no')
}
print(arg)
#arg <- c( "/home/wmj/02_Metagenome/16s/20191012_CYM/Analysis/10_LEfSe/Total/Total.kegg.L3.res","/home/wmj/02_Metagenome/16s/20191012_CYM/Analysis/17_kegg_L3/Total" , "/home/wmj/02_Metagenome/16s/20191012_CYM/bin/Total-color.txt" )
setwd(arg[2])
library(RColorBrewer)
data <- read.table(arg[1],row.names=1,sep="\t")

grp <- read.table(arg[3],header=T,sep="\t", comment.char = "")
rownames(grp)<-as.character(grp[,1])

gg <- unique(as.character(grp[,2]))
cc <-c()
for(i in 1:length(gg)){
	g1 <- grep(paste0("^",gg[i],"$"),as.character(data[,2]))
	cc <-c(cc,g1)

}

data <-data[cc,]

col <-unique(as.character(grp[,3]))
names(col)=gg
dd <- data[,2:3]
ord <- order(as.character(dd[,1]),as.numeric(dd[,2]))
dd <-dd[ord,]

gg <-as.character(unique(dd[,1]))
cc <-c()
for(i in 1:length(gg)){
	n1 <- grep(paste0("^",gg[i],"$"),as.character(dd[,1]))
	cc <- c(cc,rep(col[i],length(n1)))


}

dd <- cbind(dd,cc)

f1 <- basename(arg[1])
f2 <- strsplit(f1,".",fixed=T)[[1]][1]
h0=length(dd[,2])
if(h0>15){
	hh=ceiling(h0/4)
} else {
	hh=h0
}
hh = 5
pdf(paste0(f1,".kegg.L3.pdf"),h=hh,w=hh*1.5)
par(mar=c(3,20,4,3),cex.axis=1.5)
pos <- barplot(as.numeric(dd[,2]),horiz=TRUE,col=as.character(dd[,3]),width=0.5)
if(length(gg) > 4){
	legend("top",legend=gg,fill=unique(cc),bty="n",cex=1.2,inset=-0.1,xpd=T,ncol=4)
} else {
	legend("top",legend=gg,fill=unique(cc),bty="n",horiz=T,cex=1.2,inset=-0.1,xpd=T)
}
mtext(text=rownames(dd),side=2,adj=1,at=pos,las=1,cex=1.2,xpd=T)
dev.off()

png(paste0(f1,".kegg.L3.png"),h=hh*480/7,w=hh*1.5*480/7, bg="transparent")
par(mar=c(3,20,4,3),cex.axis=1.5)
pos <- barplot(as.numeric(dd[,2]),horiz=TRUE,col=as.character(dd[,3]),width=0.5)
if(length(gg) > 4){
	legend("top",legend=gg,fill=unique(cc),bty="n",cex=1.2,inset=-0.1,xpd=T,ncol=4)
} else {
	legend("top",legend=gg,fill=unique(cc),bty="n",horiz=T,cex=1.2,inset=-0.1,xpd=T)
}
mtext(text=rownames(dd),side=2,adj=1,at=pos,las=1,cex=1.2,xpd=T)
dev.off()


