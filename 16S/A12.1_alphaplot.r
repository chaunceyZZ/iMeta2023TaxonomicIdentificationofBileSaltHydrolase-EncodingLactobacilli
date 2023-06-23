arg <- commandArgs(T)
if(length(arg) < 1){
	cat("Argument: Out_Dir \n")
	quit('no')
}
print(arg)
#arg <- c( "/home/wmj/02_Metagenome/16s/20191012_CYM/Analysis/7_OTUs_alpha/Total/Total.final.otu.simpson.txt","/home/wmj/02_Metagenome/16s/20191012_CYM/bin/Total-color.txt"   , "/home/wmj/02_Metagenome/16s/20191012_CYM/Analysis/16_alphaRarefactionPlot/Total"  )

data <- read.table(arg[2],header=F,na.strings="NA",fill=TRUE,check.names = FALSE) ##fill=T,含有空的时候
#grp <-read.table(arg[3],header=T,row.names=1,check.names = FALSE)
grp <-read.table(arg[3],header=T,check.names = FALSE,row.names=1, comment.char = "")





SamID <- strsplit(basename(arg[2]),'.',fixed=T)[[1]][1] ####### NOTE
grp.nn <- seq(3,nrow(data)-1,by=3)

#########绘图
#col.err=c("#DC143C64","#4169E164","#2E8B5764","#9932CC64","#FF8C0064","#FFFF0064","#8B451364","#FF69B464","#80808064")
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=unique(as.character(grp[,2]))
names(col.line)=unique(as.character(grp[,1]))
col.err=col.line

### 自定义标准误
plot.error <- function(x, y, sd, len = 1, col = "black",lwd=1,lty=1) {
	len <- len * 0.05 
	arrows(x0 = x, y0 = y, x1 = x, y1 = y - sd, col = col, angle = 90, length = len,lty=lty,lwd=lwd) 
	arrows(x0 = x, y0 = y, x1 = x, y1 = y + sd, col = col, angle = 90, length = len,lty=lty,lwd=lwd) 
}

### 坐标轴
x <- data[1,-1]
y <- data[grp.nn+1,-1]
grp.name <- data[grp.nn,1]
err <-data[grp.nn+2,-1]
x.max <- data[2,2]
y1 <- unlist(y)
err1 <- unlist(err)
n.max <- max(y1[!is.na(y1)])
y.max <- n.max +err1[which(y1 == n.max)]
y.max <- ceiling(y.max)*1.5
tt1 <- strsplit(SamID,"Group")[[1]][1]

f1 <- strsplit(basename(arg[3]),"-",fixed=T)[[1]]
f2 <- f1[1]
dir <- paste(arg[1],"/",f2,sep="")
if ( !file.exists(dir))	dir.create(dir) #创建名为tmp的文件夹 	
setwd(dir)	

### plot

ggk <- grp
ggk1 <- 1
un.k <- as.character(unique(grp[,1]))
leg <- c()
#ylab <-paste("Rarefaction Measure:",tt1,sep="")
ylab=tt1
#main <- paste(tt1,":Group(",un.gg[k],")",sep="")
pdf(paste(f1[1],".",SamID,".rarefaction.pdf",sep=""))
par(mar = c(12,14,7,6),mex=0.5,cex=1,font=1,xaxs = "i", yaxs = "i",adj=0.5,mgp=c(9,1.5,0),lty=1,lwd=2,cex.lab=2.5,cex.axis=2,tck=0.02)
plot(x=0,y=0,col="white",xlim=c(0,x.max),ylim=c(0,y.max),xlab="",ylab="",las=1) #las=1:y轴刻度水平
for (l in 1:length(un.k)){ 
	nn <-as.character(unique( ggk[which(ggk[,ggk1]==un.k[l]),1]))
	n1 <- grep(paste("^",as.character(nn),"$",sep=""),as.character(grp.name))
	#leg <- c(leg,as.character(grp.name[n1]))
	x1 <- x
	y1 <- y[n1,]
	err1 <- err[n1,]
	tt <- which(is.na(y1))
	if (length(tt) >0){	### 标准误
		#lines(x=as.numeric(x1[-tt]),y=as.numeric(y1[-tt]),col=col[l])
		ee <- as.numeric(err1[-tt])
		nne <- seq(0,length(ee),length.out=10)
		plot.error(x=as.numeric(x1[-tt])[nne],y=as.numeric(y1[-tt])[nne],sd=ee[nne],col=col.err[l],lwd=1,lty=1)	
	} else {
		#lines(x=as.numeric(x1),y=as.numeric(y1),col=col[l])
		ee <- as.numeric(err1)
		nne <- seq(0,length(ee),length.out=10)
		plot.error(x=as.numeric(x1)[nne],y=as.numeric(y1)[nne],sd=ee[nne],col=col.err[l],lwd=1,lty=1)
	}
				
}
for (l in 1:length(un.k)){ ### 画曲线
	nn <-unique( ggk[which(ggk[,ggk1]==un.k[l]),1])
	n1 <- grep(paste("^",as.character(nn),"$",sep=""),as.character(grp.name))
	leg <- c(leg,as.character(grp.name[n1]))
	x1 <- x
	y1 <- y[n1,]
	err1 <- err[n1,]
	tt <- which(is.na(y1))
	if (length(tt) >0){	### 曲线
		lines(x=as.numeric(x1[-tt]),y=as.numeric(y1[-tt]),col=col.line[l])
		#plot.error(x=as.numeric(x1[-tt]),y=as.numeric(y1[-tt]),sd=as.numeric(err1[,-tt]),col=col[l],lwd=1,lty=1)	
	} else {
		lines(x=as.numeric(x1),y=as.numeric(y1),col=col.line[l])
		#plot.error(x=as.numeric(x1),y=as.numeric(y1),sd=as.numeric(err1),col=col[l],lwd=0.5,lty=1)
	}
				
}
mtext(text="Sequences Per Sample",side=1,line=7,adj=0.5,cex=2.5)
mtext(text=ylab,side=2,line=9,adj=0.5,cex=2.5)
legend("topleft",legend=leg,ncol=1,cex=2,bty="n",col=col.line[1:l],lty=1)
dev.off()
		
png(paste(f1[1],".",SamID,".rarefaction.png",sep=""), bg="transparent")
leg=c()
par(mar = c(12,14,7,6),mex=0.5,cex=1,font=1,xaxs = "i", yaxs = "i",adj=0.5,mgp=c(9,1.5,0),lty=1,lwd=2,cex.lab=2.5,cex.axis=2,tck=0.02)
plot(x=0,y=0,col="white",xlim=c(0,x.max),ylim=c(0,y.max),xlab="",ylab="",las=1) #las=1:y轴刻度水平
for (l in 1:length(un.k)){ 
	nn <-as.character(unique( ggk[which(ggk[,ggk1]==un.k[l]),1]))
	n1 <- grep(paste("^",as.character(nn),"$",sep=""),as.character(grp.name))
	#leg <- c(leg,as.character(grp.name[n1]))
	x1 <- x
	y1 <- y[n1,]
	err1 <- err[n1,]
	tt <- which(is.na(y1))
	if (length(tt) >0){	### 标准误
		#lines(x=as.numeric(x1[-tt]),y=as.numeric(y1[-tt]),col=col[l])
		ee <- as.numeric(err1[-tt])
		nne <- seq(0,length(ee),length.out=10)
		plot.error(x=as.numeric(x1[-tt])[nne],y=as.numeric(y1[-tt])[nne],sd=ee[nne],col=col.err[l],lwd=1,lty=1)	
	} else {
		#lines(x=as.numeric(x1),y=as.numeric(y1),col=col[l])
		ee <- as.numeric(err1)
		nne <- seq(0,length(ee),length.out=10)
		plot.error(x=as.numeric(x1)[nne],y=as.numeric(y1)[nne],sd=ee[nne],col=col.err[l],lwd=1,lty=1)
	}
				
}
for (l in 1:length(un.k)){ ### 画曲线
	nn <-unique( ggk[which(ggk[,ggk1]==un.k[l]),1])
	n1 <- grep(paste("^",as.character(nn),"$",sep=""),as.character(grp.name))
	leg <- c(leg,as.character(grp.name[n1]))
	x1 <- x
	y1 <- y[n1,]
	err1 <- err[n1,]
	tt <- which(is.na(y1))
	if (length(tt) >0){	### 曲线
		lines(x=as.numeric(x1[-tt]),y=as.numeric(y1[-tt]),col=col.line[l])
		#plot.error(x=as.numeric(x1[-tt]),y=as.numeric(y1[-tt]),sd=as.numeric(err1[,-tt]),col=col[l],lwd=1,lty=1)	
	} else {
		lines(x=as.numeric(x1),y=as.numeric(y1),col=col.line[l])
		#plot.error(x=as.numeric(x1),y=as.numeric(y1),sd=as.numeric(err1),col=col[l],lwd=0.5,lty=1)
	}
				
}
mtext(text="Sequences Per Sample",side=1,line=7,adj=0.5,cex=2.5)
mtext(text=ylab,side=2,line=9,adj=0.5,cex=2.5)
legend("topleft",legend=leg,ncol=1,cex=2,bty="n",col=col.line[1:l],lty=1)
dev.off()
		
		




 
