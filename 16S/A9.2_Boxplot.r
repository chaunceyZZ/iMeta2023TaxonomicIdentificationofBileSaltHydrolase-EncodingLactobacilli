arg <- commandArgs(T)
if(length(arg) != 3){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}
print(arg)
#arg <-c("/home/wmj/02_Metagenome/16s/20191017_LY/Analysis/9_PICRUSt/Total/Total.kegg.L3.abundance.txt", "/home/wmj/02_Metagenome/16s/20191017_LY/Analysis/10_LEfSe/Total" , "/home/wmj/02_Metagenome/16s/20191017_LY/Analysis/9_PICRUSt/Total/Total.catagrized.kegg.L2.tab")
library(RColorBrewer)
setwd(arg[2])

ff=unlist(strsplit(basename(arg[1]),'\\.'))[1]

abundance=read.table(arg[1],sep='\t',header=T,row.names=1)
gr=read.table(arg[3],sep='\t',header=T,row.names=1,comment.char = "")
grp=matrix(NA,nrow=nrow(gr),ncol=2)
colnames(grp)=c('level','KEGG')
for(i in 1:nrow(gr)){
	ii=as.character(gr[i,ncol(gr)])
	ii1=unlist(strsplit(ii,'; '))
	grp[i,1]=ii1[1]
	grp[i,2]=ii1[2]
}
rownames(grp)=grp[,2]


data=cbind(abundance,grp[rownames(abundance),1])
or=order(data[,ncol(data)],rownames(data),decreasing=T)
dd=data[or,]
colnames(dd)=c(colnames(abundance),'gr')

#col=brewer.pal(9,"Pastel1")
col=brewer.pal(12,"Set3")[c(1,3:7,9:12)]
pdf(paste0(ff,'.kegg.boxplot.pdf'),h=8,w=10)
lay_text <- matrix(c(1:2,0,3),ncol=2,nrow=2,byrow=T)
layout(mat=lay_text,w=c(1.5,1),h=c(12,1))	

############# legend
par(mar=c(2,0,2,0))
plot(0,0,xlab="",ylab="" ,axes=F ,xlim=c(0,10),ylim=c(0,nrow(dd)),col="white")
tt=unique(as.character(dd[,'gr']))
yy0=0
for(i in 1:length(tt)){
	i1=as.character(dd[,'gr']) %in% tt[i]
	d1=dd[i1,]
	yy0 <-yy0 -0.5
	yy1 <-yy0+nrow(d1) +0.5
	rect(xleft=0, ybottom=yy0, xright=10,ytop=yy1,col=col[i],border=col[i],xpd=T)
	text(x=10,y=seq(yy0+1,yy1-0.5,length.out=nrow(d1))-0.5,labels=rownames(d1),adj=1)
#	text(x=0,y=yy1-1,labels=tt[i],adj=0)
	i2=gsub('\\s','\\\n',tt[i])
	text(x=0,y=yy1-1,labels=tt[i],adj=0,font=2)
	yy0=yy1
}

########boxplot
yy0=0
par(mar=c(2,0,2,2))
plot(0,0,xlab="",ylab="",axes=F,xlim=c(min(abundance),max(abundance)),ylim=c(0,nrow(dd)),col="white") #axes=F
for(i in 1:length(tt)){
	i1=as.character(dd[,'gr']) %in% tt[i]
	d1=dd[i1,-ncol(dd)]
	d2=t(d1)
	yy0 <-yy0 -0.5
	yy1 <-yy0+nrow(d1) +0.5
	yy=seq(yy0+1,yy1-0.5,length.out=nrow(d1))-0.5
	boxplot(horizontal=T,col=col[i],d2,add=T,axes=F,at=yy,cex=0.5,pch=16,border=col[i])
	yy0=yy1
}
####### legend
par(mar=c(0,0,0,2))
plot(0,0,xlab="",ylab="",axes=F,xlim=c(min(abundance),max(abundance)),ylim=c(0,2),col="white") #axes=F
ll=seq(min(abundance),max(abundance),length.out=5)
axis(side=1,line=-6,at=round(ll,2),labels=round(ll,2))
mtext(side=1,line=-3.5,xpd=T,text="Relative Abundance")
dev.off()

##########
###########
###########
png(paste0(ff,'.kegg.boxplot.png'),h=8*480/7,w=10*480/7)
lay_text <- matrix(c(1:2,0,3),ncol=2,nrow=2,byrow=T)
layout(mat=lay_text,w=c(1.5,1),h=c(12,1))	

############# legend
par(mar=c(2,0,2,0))
plot(0,0,xlab="",ylab="" ,axes=F ,xlim=c(0,10),ylim=c(0,nrow(dd)),col="white")
tt=unique(as.character(dd[,'gr']))
yy0=0
for(i in 1:length(tt)){
	i1=as.character(dd[,'gr']) %in% tt[i]
	d1=dd[i1,]
	yy0 <-yy0 -0.5
	yy1 <-yy0+nrow(d1) +0.5
	rect(xleft=0, ybottom=yy0, xright=10,ytop=yy1,col=col[i],border=col[i],xpd=T)
	text(x=10,y=seq(yy0+1,yy1-0.5,length.out=nrow(d1))-0.5,labels=rownames(d1),adj=1)
#	text(x=0,y=yy1-1,labels=tt[i],adj=0)
	i2=gsub('\\s','\\\n',tt[i])
	text(x=0,y=yy1-1,labels=tt[i],adj=0,font=2)
	yy0=yy1
}

########boxplot
yy0=0
par(mar=c(2,0,2,2))
plot(0,0,xlab="",ylab="",axes=F,xlim=c(min(abundance),max(abundance)),ylim=c(0,nrow(dd)),col="white") #axes=F
for(i in 1:length(tt)){
	i1=as.character(dd[,'gr']) %in% tt[i]
	d1=dd[i1,-ncol(dd)]
	d2=t(d1)
	yy0 <-yy0 -0.5
	yy1 <-yy0+nrow(d1) +0.5
	yy=seq(yy0+1,yy1-0.5,length.out=nrow(d1))-0.5
	boxplot(horizontal=T,col=col[i],d2,add=T,axes=F,at=yy,cex=0.5,pch=16,border=col[i])
	yy0=yy1
}
####### legend
par(mar=c(0,0,0,2))
plot(0,0,xlab="",ylab="",axes=F,xlim=c(min(abundance),max(abundance)),ylim=c(0,2),col="white") #axes=F
ll=seq(min(abundance),max(abundance),length.out=5)
axis(side=1,line=-6,at=round(ll,2),labels=round(ll,2))
mtext(side=1,line=-3.5,xpd=T,text="Relative Abundance")
dev.off()


