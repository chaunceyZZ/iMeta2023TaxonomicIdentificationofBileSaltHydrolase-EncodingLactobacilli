arg <- commandArgs(T)
if(length(arg) > 5){
        cat("Argument: Out_Dir Data_File\n")
        quit('no')
}
print(arg)
#arg <- c("/home/wmj/02_Metagenome/16s/HTQ_190723/20190904/Analysis/1_Abundance"  , "/home/wmj/02_Metagenome/16s/HTQ_190723/20190904/Analysis/1_Abundance/Total.PlyGen.abundance.txt", "Total.Phylum.abundance.txt"  ,"/home/wmj/02_Metagenome/16s/HTQ_190723/20190904/bin/Total-color.txt")

setwd(arg[1])
library(RColorBrewer)

grpInfo <- read.table(arg[4],header=T,sep="\t",check.names = FALSE, comment.char = "")
rownames(grpInfo) <- grpInfo[,1]


####	Genus
PlyGen <- read.table(arg[2],header=T,row.names=1,sep="\t",check.names = FALSE)
#####	Phylum

Phylum <- read.table(arg[3],header=T,sep="\t",check.names = FALSE)

i1=intersect(rownames(grpInfo),colnames(PlyGen))
ii <- intersect(i1,colnames(Phylum))
if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),ii)
	print(ii1)
}
grpInfo=grpInfo[ii,]
### 去除全NA的列
g.na <- apply(PlyGen,2,sum)
g.nn <- which(is.na(g.na))
if(length(g.nn)!=0) PlyGen <- PlyGen[,-g.nn]

PlyGen <- PlyGen[,ii]
### 拆分
row.ply <- rownames(PlyGen)
strain.ply <- matrix(NA,nrow=nrow(PlyGen),ncol=2)
colnames(strain.ply) <- c("Ply","Gen")
for (i in 1:length(row.ply)){
	tt <- unlist(strsplit(row.ply[i],split=" "))
	strain.ply[i,] <-tt
}
rownames(PlyGen) <- strain.ply[,2]

### 	Top.5
top.ply =10
ply <- cbind(PlyGen,rowSums(PlyGen));
strain.ply <-strain.ply[order(ply[,ncol(ply)],decreasing=T),];
ply <- ply[order(ply[,ncol(ply)],decreasing=T),];
if(nrow(ply) > top.ply){	
	strain.ply <- strain.ply[-c((top.ply+1):nrow(ply)),];
	ply <- ply[-c((top.ply+1):nrow(ply)),];
}
ply <- ply[,-ncol(ply)];
ord <- order(strain.ply[,1],strain.ply[,2],decreasing=T)
strain.ply<-strain.ply[ord,]
ply <- ply[ord,]


## 计算比例
sum.ply <- apply(ply,2,sum) 
sum1 <- matrix( rep(sum.ply,nrow(ply)) , nrow=nrow(ply),byrow=T)
Ply <- as.matrix(ply/sum1) 

### color
#newpalette<-colorRampPalette(brewer.pal(9,"Blues"))(10)
#col=brewer.pal(12,"Paired")
col=c("Oranges","Greens","Reds","Blues","Purples","PuOr","PuBuGn","PuBu","OrRd")
ply.color<-strain.ply[,1]
names(ply.color) <-ply.color
un.ply <-unique(strain.ply[,1])
for(i in 1:length(un.ply)){
	nn <- grep(un.ply[i],strain.ply[,1])
	ply.color[nn] <- brewer.pal(9,col[i])[2:(1+length(nn))]

}


#####################
#####	Phylum

#Phylum <- read.table(arg[3],header=T,sep="\t",check.names = FALSE)

### 合并同一个门的丰度
row.phy <- Phylum[,1]

un.phy <- unique(row.phy)
phy <- matrix(NA,ncol=ncol(Phylum[,-1]),nrow=length(un.phy))
colnames(phy) <- colnames(Phylum[,-1])
rownames(phy) <- un.phy
for (i in 1:length(un.phy)){
	# nn <- grep(paste("^",as.character(un.phy[i]),"$",sep=""),as.character(Phylum[,1]))
        nn <- which(Phylum$ID %in% un.phy[i])
	if(length(nn) >1){
		phy[i,] <- apply(Phylum[nn,-1],2,sum)
	} else {
		phy[i,] <- as.matrix(Phylum[nn,-1],byrow=T,nrow=1)
	}

}
Phylum <- phy
Phylum <- Phylum[,ii]
### 取总丰度Top 5
top.phy =5
phy <- cbind(Phylum,rowSums(Phylum));
phy <- phy[order(phy[,ncol(phy)],decreasing=T),];
if(nrow(phy) > top.phy){
	phy <- phy[-c((top.phy+1):nrow(phy)),];
}
phy <- phy[,-ncol(phy)];
phy <- phy[order(rownames(phy),decreasing=T),];

## 计算比例
sum.phy <- apply(phy,2,sum) 
sum1 <- matrix( rep(sum.phy,nrow(phy)) , nrow=nrow(phy),byrow=T)
Phy <- as.matrix(phy/sum1)  


####	Group
grp.id <- as.character(grpInfo[,1])
grp.grp <- as.character(grpInfo[,2])
grp <- as.character(unique(grp.grp))
print(grp)
ax.grp <- c()
plot.grp <- c()
plot.nn <-c()
plot.at.lab <- c()
ord <- c()
for (i in 1:length(grp)){
#	nn <- grep(grp[i],grp.grp)
	nn <- grep(paste("^",grp[i],"$",sep=""),grp.grp)
	id <- grp.id[nn]
	plot.nn <- c(plot.nn,id)
	plot.at.lab <- c(plot.at.lab,c(1:length(id)))
	ax.grp <-c(ax.grp,unique(as.character(grpInfo[nn,2])))
	ord <- c(ord,nn)
	if (i==1){
		plot.grp <- c(rep(0,length(id)),1)
	} else {
		plot.grp <- c(plot.grp,rep(0,length(id)-1),1)
	}
}

proportion <- c()
for (i in 1:length(grp)){
#	nn <- grep(grp[i],grp.grp)
	nn <- grep(paste("^",grp[i],"$",sep=""),grp.grp)
#	dd <- Phy[,nn]    ## Normalization
	dd <- phy[,nn]    # Non_Nor
	n1 <- sum(as.numeric(dd["Bacteroidetes",]))/sum(as.numeric(dd["Firmicutes",]))
	n2 <- signif(n1,2)
	proportion <- c(proportion,n2 )
}
names(proportion) <- grp




SamID <- strsplit(basename(arg[2]),'.',fixed=T)[[1]][1] ####### NOTE




options(scipen=3) ## 取消科学计数法

#### Barplot


ff <- paste(SamID,".abundance.barplot.pdf",sep="")
pdf(ff,height=7.5,width=10)
layout(mat=matrix(c(1,2,3,4),nrow=2,byrow=T),heights=c(1,1),widths=c(4,1))
## Genus
par(mar = c(3,3,8,0),mex=0.5,cex=1,lwd=0.1,font=1,xaxs = "i", yaxs = "i",adj=0.5,mgp=c(2,0,0))
m.gen <- paste("TOP.10 Genus",sep="")
gen.p <- barplot(as.matrix(Ply),beside=F,axes=F,axisnames=F,width=1,space=plot.grp,col=ply.color,border="black",main=m.gen)
a1 <- gen.p[-length(gen.p)]
a2 <- length(a1)/length(grp)
a3 <- seq(from=median(a1[1:a2]),by=a2+1,length.out=length(grp))
axis(side=1,at=a3,labels=ax.grp,tick=F)
axis(side=3,at=a1,labels=plot.at.lab,tick=F)
#axis(side=2,at=1,labels="A",tick=F,las=1,line=1)
##	legend
par(mar = c(3,0,8,1),mex=0.5,cex=1,lwd=1,font=1,xaxs = "i", yaxs = "i",adj=0.5)
plot(0,0,axes=F,col="white",xlab="",ylab="")
legend("left", legend=rev(rownames(Ply)),bty="n",ncol=1,cex=0.8,fill=rev(ply.color))

### Phylum
par(mar = c(5,3,8,0),mex=0.5,cex=1,lwd=0.1,font=1,xaxs = "i", yaxs = "i",adj=0.5,mgp=c(2,0,0))
col=brewer.pal(12,"Paired")

col.phy <- col[1:nrow(Phy)]
m.phy <- paste("TOP.5 Phylum",sep="")
phy.p <- barplot(as.matrix(Phy),beside=F,axes=F,axisnames=F,width=1,space=plot.grp,col=col.phy,border="black",main=m.phy)
a1 <- phy.p[-length(phy.p)]
a2 <- length(a1)/length(grp)
a3 <- seq(from=median(a1[1:a2]),by=a2+1,length.out=length(grp))
axis(side=1,at=a3,labels=ax.grp,tick=F)
axis(side=3,at=a1,labels=plot.at.lab,tick=F)
#axis(side=2,at=1,labels="B",tick=F,las=1,line=1)
mtext(side=1,text=paste("B : F = ",proportion,sep=""),line=2,at=a3)
##	legend
par(mar = c(3,0,8,1),mex=0.5,cex=1,lwd=1,font=1,xaxs = "i", yaxs = "i",adj=0.5)
plot(0,0,axes=F,col="white",xlab="",ylab="")
legend("left", legend=rev(rownames(Phy)),bty="n",ncol=1,cex=0.8,fill=rev(col.phy))

dev.off()



ff0 <- paste(SamID,".abundance.barplot.png",sep="")
png(ff0,height=7.5*480/7,width=10*480/7, bg="transparent")
layout(mat=matrix(c(1,2,3,4),nrow=2,byrow=T),heights=c(1,1),widths=c(4,1))
## Genus
par(mar = c(3,3,8,0),mex=0.5,cex=1,lwd=0.1,font=1,xaxs = "i", yaxs = "i",adj=0.5,mgp=c(2,0,0))
m.gen <- paste("TOP.10 Genus",sep="")
gen.p <- barplot(as.matrix(Ply),beside=F,axes=F,axisnames=F,width=1,space=plot.grp,col=ply.color,border="black",main=m.gen)
a1 <- gen.p[-length(gen.p)]
a2 <- length(a1)/length(grp)
a3 <- seq(from=median(a1[1:a2]),by=a2+1,length.out=length(grp))
axis(side=1,at=a3,labels=ax.grp,tick=F)
axis(side=3,at=a1,labels=plot.at.lab,tick=F)
#axis(side=2,at=1,labels="A",tick=F,las=1,line=1)
##	legend
par(mar = c(3,0,8,1),mex=0.5,cex=1,lwd=1,font=1,xaxs = "i", yaxs = "i",adj=0.5)
plot(0,0,axes=F,col="white",xlab="",ylab="")
legend("left", legend=rev(rownames(Ply)),bty="n",ncol=1,cex=0.8,fill=rev(ply.color))

### Phylum
par(mar = c(5,3,8,0),mex=0.5,cex=1,lwd=0.1,font=1,xaxs = "i", yaxs = "i",adj=0.5,mgp=c(2,0,0))
col=brewer.pal(12,"Paired")

col.phy <- col[1:nrow(Phy)]
m.phy <- paste("TOP.5 Phylum",sep="")
phy.p <- barplot(as.matrix(Phy),beside=F,axes=F,axisnames=F,width=1,space=plot.grp,col=col.phy,border="black",main=m.phy)
a1 <- phy.p[-length(phy.p)]
a2 <- length(a1)/length(grp)
a3 <- seq(from=median(a1[1:a2]),by=a2+1,length.out=length(grp))
axis(side=1,at=a3,labels=ax.grp,tick=F)
axis(side=3,at=a1,labels=plot.at.lab,tick=F)
#axis(side=2,at=1,labels="B",tick=F,las=1,line=1)
mtext(side=1,text=paste("B : F = ",proportion,sep=""),line=2,at=a3)
##	legend
par(mar = c(3,0,8,1),mex=0.5,cex=1,lwd=1,font=1,xaxs = "i", yaxs = "i",adj=0.5)
plot(0,0,axes=F,col="white",xlab="",ylab="")
legend("left", legend=rev(rownames(Phy)),bty="n",ncol=1,cex=0.8,fill=rev(col.phy))

dev.off()

################

#### Barplot


ff <- paste(SamID,".abundance.barplot.tmp.pdf",sep="")
pdf(ff,height=7.5,width=10)
layout(mat=matrix(c(1,2,3,4),nrow=2,byrow=T),heights=c(1,1),widths=c(4,1))
## Genus
par(mar = c(3,3,8,0),mex=0.5,cex=1,lwd=0.1,font=1,xaxs = "i", yaxs = "i",adj=0.5,mgp=c(2,0,0))
m.gen <- paste("TOP.10 Genus",sep="")
t1=as.matrix(Ply)
tt=matrix(NA,nrow=nrow(t1),ncol=length(grp))
rownames(tt)=rownames(t1)
colnames(tt)=grp
for(j in grp){
	g1=as.character(grpInfo[,2]) %in% j
	g2=as.character(grpInfo[g1,1])
	tt[,j]=apply(t1[,g2],1,mean)
}
gen.p <- barplot(tt,beside=F,axes=F,axisnames=F,width=1,col=ply.color,border="black",main=m.gen)
axis(side=1,at=gen.p,labels=ax.grp,tick=F)
#axis(side=2,at=1,labels="A",tick=F,las=1,line=1)
##	legend
par(mar = c(3,0,8,1),mex=0.5,cex=1,lwd=1,font=1,xaxs = "i", yaxs = "i",adj=0.5)
plot(0,0,axes=F,col="white",xlab="",ylab="")
legend("left", legend=rev(rownames(Ply)),bty="n",ncol=1,cex=0.8,fill=rev(ply.color))

### Phylum
par(mar = c(5,3,8,0),mex=0.5,cex=1,lwd=0.1,font=1,xaxs = "i", yaxs = "i",adj=0.5,mgp=c(2,0,0))
col=brewer.pal(12,"Paired")

col.phy <- col[1:nrow(Phy)]
m.phy <- paste("TOP.5 Phylum",sep="")

t1=as.matrix(Phy)
tt=matrix(NA,nrow=nrow(t1),ncol=length(grp))
rownames(tt)=rownames(t1)
colnames(tt)=grp
for(j in grp){
	g1=as.character(grpInfo[,2]) %in% j
	g2=as.character(grpInfo[g1,1])
	tt[,j]=apply(t1[,g2],1,mean)
}


phy.p <- barplot(tt,beside=F,axes=F,axisnames=F,width=1,col=col.phy,border="black",main=m.phy)
axis(side=1,at=phy.p,labels=ax.grp,tick=F)
#axis(side=2,at=1,labels="B",tick=F,las=1,line=1)
mtext(side=1,text=paste("B : F = ",proportion,sep=""),line=2,at=a3)
##	legend
par(mar = c(3,0,8,1),mex=0.5,cex=1,lwd=1,font=1,xaxs = "i", yaxs = "i",adj=0.5)
plot(0,0,axes=F,col="white",xlab="",ylab="")
legend("left", legend=rev(rownames(Phy)),bty="n",ncol=1,cex=0.8,fill=rev(col.phy))

dev.off()

################





