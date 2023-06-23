

arg <- commandArgs(T)
if(length(arg) != 5){
        cat("Argument: Data_File Group_File Out_Dir Group_column Color\n")
        quit('no')
}
print(arg)
#arg <- c( "/home/Project/Metabolome/HFQ/HFQ_LOESS/HILIC/0.0_raw_data/FAP.data.txt","/home/Project/Metabolome/HFQ/HFQ_LOESS/HILIC/bin/FAP-grouping.info", "/home/Project/Metabolome/HFQ/HFQ_LOESS/HILIC/3.2_OPLS_DA", "2" ,"/home/Project/Metabolome/HFQ/HFQ_LOESS/HILIC/bin/FAP-color.txt" )
library("ropls")
library("ggplot2")
library("car")
library("rgl");
library("scatterplot3d");
library("plot3D")

dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);

#dataFile <- "FT.data.txt";
#groupFile <- "FT-grouping.info";
#outDir <- "./";
#grp_col <- 2;

grpInfo <- read.table(groupFile,header=T);
FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- FileName[1];
Data <- read.table(dataFile,header=T,sep="\t");
nn <- as.character(Data[,1])
dd <- as.matrix(Data[,-1])
rownames(dd) <-nn
Data <-dd
Data.tm <- t(as.matrix(Data)); 
groupname <- c();

for(i in 1:length(colnames(Data)))
{
	groupname <- append(groupname,as.character(grpInfo[grep(paste("^",colnames(Data)[i],"$",sep=""),as.character(grpInfo[,1])),grp_col]));
}

Groups <- unique(groupname)

if(length(Groups) == 2){
	Data.oplsda <- try(opls(Data.tm ,groupname, predI =1, orthoI = NA,plotL = FALSE, printL =FALSE),silent=TRUE);
	if ('try-error' %in% class(Data.oplsda))
	{
		Data.oplsda <- opls(Data.tm ,groupname, predI =1, orthoI = 2,plotL = FALSE, printL =FALSE);
	} 
	plot(Data.oplsda,typeVc="summary")
	plot(Data.oplsda,typeVc="x-score",parLabVc=groupname)
	VIP <- getVipVn(Data.oplsda);
	write.table(as.matrix(VIP),file=paste(outDir,"/",SamID,".All.OPLSDA.VIP.txt",sep=""),quote = F,sep="\t",col.names=NA,row.names=T);
	VIP <- VIP[VIP>1];
	write.table(as.matrix(VIP),file=paste(outDir,"/",SamID,".OPLSDA.VIP.txt",sep=""),quote = F,sep="\t",col.names=NA,row.names=T);
}else if(length(Groups) >2){
	Data.oplsda <- opls(Data.tm ,groupname,predI =1, plotL = FALSE, printL =FALSE);
        plot(Data.oplsda,typeVc="summary")
	plot(Data.oplsda,typeVc="x-score",parLabVc=groupname)
        VIP <- getVipVn(Data.oplsda);
	write.table(as.matrix(VIP),file=paste(outDir,"/",SamID,".All.OPLSDA.VIP.txt",sep=""),quote = F,sep="\t",col.names=NA,row.names=T);
        VIP <- VIP[VIP>1];
        write.table(as.matrix(VIP),file=paste(outDir,"/",SamID,".OPLSDA.VIP.txt",sep=""),quote = F,sep="\t",col.names=NA,row.names=T);
}

dev.off()


if(length(Groups) == 2){
	####	Points
	dd <- attributes(Data.oplsda)
	nn <- attributes(dd)
	dd.s <-dd$suppLs
	nn.s <- attributes(dd.s)
	lab <-dd$summaryDF
	p.info <- read.table(arg[5],header=T,row.names=1,comment.char = "")
	#p.info[,2] <- paste("#",p.info[,2],sep="")
	p.nn <- rownames(dd$scoreMN)
	ii <- intersect(p.nn,rownames(p.info))
	p.info <- p.info[ii,]
	#pch <- as.numeric(p.info[,2])
	x <- dd$scoreMN
	y <- dd$orthoScoreMN[,1]
	
	col <- as.character(p.info[,2])
	pdf(paste(outDir,"/",SamID,".OPLSDA.point.pdf",sep=""),height=8,width=8)
	par(mgp=c(2.5,1,0),mar=c(10,5,2,2))
	pp <-signif(dd$modelDF[1,1],2)*100
	plot(x=x,y=y,pch=rep(20,length(col)),col=col,main=paste(SamID," Scores(OPLS-DA)",sep=""),xlab=paste("t1(",pp,"%)",sep=""),ylab="to1",cex.lab=2,cex.main=2,cex.axis=2,xlim=c(-30,30),ylim=c(-40,40))

#dataEllipse(plot[grep(Groups[1],groupname),],levels = c(0.95), add=TRUE, col = brewer.pal(8,"Set2")[1], lwd = 1,plot.points=FALSE,fill=TRUE,center.cex=0.2)
	dataEllipse(x=as.numeric(x),y=as.numeric(y),groups=factor(groupname,levels=unique(groupname),labels=unique(groupname)),levels = c(0.95), add=TRUE, col =unique(col) , lwd = 1,plot.points=FALSE,fill=TRUE,center.cex=0.2,group.labels="",fill.alpha = 0)
#abline(h=0,lwd=2)
#abline(v=0,lwd=2)
	legend("topright",pch=rep(20,length(unique(col))),col=unique(col),bty="n",legend=as.character(unlist(unique(p.info[,1]))),cex=1.8,pt.cex=1.8)
#l1 <- c( "R2X","R2Y","Q2"," RMSEE"," pre"," ort")
#l2 <-as.numeric(lab[1:6])
#at= seq(min(x),max(x),length.out=6)
	l1 <- c( "R2X","R2Y","Q2")
	l2 <-as.numeric(lab[1:3])
#at= seq(min(x),max(x),length.out=3)
	at=seq(min(x),max(x),length.out=3)
	mtext(side=1,text=l1,at=at,line=4,adj=0.5,cex=2,font=2)
	mtext(side=1,text=l2,at=at,line=5.5,adj=0.5,cex=2,font=2)

	dev.off()

###	Area
	pdf(paste(outDir,"/",SamID,".OPLSDA.area.pdf",sep=""),height=8,width=8)
	par(mgp=c(2.5,1,0),mar=c(10,5,2,2))
	pp <-signif(dd$modelDF[1,1],2)*100
	plot(x=x,y=y,pch=rep(20,length(col)),col=col,main=paste(SamID," Scores(OPLS-DA)",sep=""),xlab=paste("t1(",pp,"%)",sep=""),ylab="to1",cex.lab=2,cex.main=2,cex.axis=2,xlim=c(-30,30),ylim=c(-40,40))

	
	dataEllipse(x=as.numeric(x),y=as.numeric(y),groups=factor(groupname,levels=unique(groupname),labels=unique(groupname)),levels = c(0.95), add=TRUE, col =unique(col) , lwd = 2,plot.points=FALSE,fill=TRUE,center.cex=0.2,group.labels="",fill.alpha = 0.3)
#abline(h=0,lwd=2)
#abline(v=0,lwd=2)
	legend("topright",pch=rep(20,length(unique(col))),col=unique(col),bty="n",legend=as.character(unlist(unique(p.info[,1]))),cex=1.8,pt.cex=1.8)
#l1 <- c( "R2X","R2Y","Q2"," RMSEE"," pre"," ort")
#l2 <-as.numeric(lab[1:6])
#at= seq(min(x),max(x),length.out=6)
	l1 <- c( "R2X","R2Y","Q2")
	l2 <-as.numeric(lab[1:3])
	at= seq(-30,30,length.out=3)
	mtext(side=1,text=l1,at=at,line=4,adj=0.5,cex=2,font=2)
	mtext(side=1,text=l2,at=at,line=5.5,adj=0.5,cex=2,font=2)

	dev.off()

###	3D
	z <- dd$orthoScoreMN[,2]	
	pdf(paste(outDir,"/",SamID,".OPLSDA.3D.pdf",sep=""),height=8,width=8)


#scatter3d(x = x, y = y, z = z, groups = as.factor(groupname),xlab=paste("t1(",pp,"%)",sep=""),ylab="to1",zlab="to2",surface=F,fill=T ,ellipsoid = TRUE,surface.col=unique(col),surface.alpha=0.3,grid = FALSE,axis.col=rep("black",3),point.col=col)

#rgl.viewpoint(zoom=0.7)

	scatter3D(x,y,z,colvar=NULL,col=col,cex=2,pch=20,xlab=paste("t1(",pp,"%)",sep=""),ylab="to1",zlab="to2",bty="b2",phi=20,theta=30,cex.lab=2)
	legend("topright",pch=rep(20,length(unique(col))),col=unique(col),bty="n",legend=as.character(unlist(unique(p.info[,1]))),cex=2,pt.cex=2,horiz=F,xpd=T,inset=-0.15)
#axes3d(edges=c("x--","y+-","z--"),cex=0.75)
	l1 <- c( "R2X","R2Y","Q2")
	l2 <-as.numeric(lab[1:3])
	legend( "bottom",horiz=T,xpd=T,inset=-0.12,legend=l1,bty="n",cex=2)
	legend( "bottom",horiz=T,xpd=T,inset=-0.17,legend=l2,bty="n",cex=2)

	dev.off()



}









