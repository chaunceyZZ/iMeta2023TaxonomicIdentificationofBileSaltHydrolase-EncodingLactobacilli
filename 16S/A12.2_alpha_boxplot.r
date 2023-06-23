arg <- commandArgs(T)
if(length(arg) < 1){
	cat("Argument: Out_Dir \n")
	quit('no')
}
print(arg)
#arg <- c( "/home/wmj/02_Metagenome/16s/20191012_CYM/Analysis/7_OTUs_alpha/Total/Total.final.otu.simpson.txt", "/home/wmj/02_Metagenome/16s/20191012_CYM/bin/Total-color.txt","/home/wmj/02_Metagenome/16s/20191012_CYM/Analysis/16_alphaRarefactionPlot/Total"   )
#print(arg)
dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- 2

grpInfo <- read.table(groupFile,header=T,check.names = FALSE, comment.char = "")
rownames(grpInfo) <- as.character(grpInfo[,1])

FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <-FileName[1];
MetID <- FileName[4];


Alpha <- as.matrix(data.frame(read.table(dataFile,header=T,row.names = 1,check.names = FALSE)));

ii <- intersect(rownames(grpInfo),rownames(Alpha))
if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),ii)
	print(ii1)
}




Alpha <- as.matrix(Alpha[ii,])
grpInfo <- grpInfo[ii,]


groupname <- c();
for(i in 1:length(rownames(Alpha)))
{
	groupname <- append(groupname,as.character(grpInfo[grep(paste("^",rownames(Alpha)[i],"$",sep=""),grpInfo[,1]),grp_col]));
}

rownames(Alpha) <- groupname;
#Groups <- as.character(levels(as.factor(grpInfo[,grp_col])))
Groups <- as.character(unique(grpInfo[,grp_col]))

#gg <- read.table(arg[5],header=T,check.names = FALSE)
#rownames(gg) <- as.character(gg[,1])
#ii <- intersect(rownames(grpInfo),rownames(gg))
#gg <- gg[ii,]
gg <-grpInfo
Groups.1 <-as.character(unique(gg[,2])) 

Alpha.Group <- c();
Alpha.Group.mean  <- c();
Alpha.Group.sd  <- c();
for(i in 1:length(Groups))
{
	Alpha.Group[[i]] <- as.matrix(Alpha[grep(paste("^",Groups[i],"$",sep=""),rownames(Alpha)),]);
	Alpha.Group.mean[[i]] <- mean (Alpha.Group[[i]]);
	Alpha.Group.sd[[i]] <- sd (Alpha.Group[[i]]);
}

Alpha.GG.Utest <- c();
Alpha.GG.Asterisk <- c();
for(m in 1:(length(Groups)-1))
{
	for(n in (m+1):length(Groups))
	{
		GG=paste(Groups[m],"-",Groups[n],sep="");
		Alpha.GG.Utest[[GG]]=wilcox.test(Alpha.Group[[m]],Alpha.Group[[n]]);
		if(Alpha.GG.Utest[[GG]]$p.value < 0.001){
			Alpha.GG.Asterisk[[GG]] <- '***' ;
		}
		else if(Alpha.GG.Utest[[GG]]$p.value < 0.01){
                        Alpha.GG.Asterisk[[GG]] <- ' **' ;
		}
		else if(Alpha.GG.Utest[[GG]]$p.value < 0.05){
                        Alpha.GG.Asterisk[[GG]] <- '  *' ;
		}
		else{
			Alpha.GG.Asterisk[[GG]] <- '   ' ;
		}
	}
}



pdf(paste(outDir,"/",SamID,".",MetID,".final.otu.barplot.pdf",sep=""),width=4,height=4)
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=unique(as.character(grpInfo[,3]))
names(col.line)=unique(as.character(grpInfo[,2]))

y.lim <- max(Alpha.Group.mean+Alpha.Group.sd)*1.3;
err.wid <- 0.1
sig.tick <- 0.03
gap <- 0.008
step <- 0.03
#par(mar = c(12,14,7,6),mex=0.5,cex=1,font=1,xaxs = "i", yaxs = "i",adj=0.5,mgp=c(9,1.5,0),lty=1,lwd=2,cex.lab=2.5,cex.axis=2,tck=0.02)
par(mgp=c(4,1,0),tck=-0.02,mar=c(2,6,2,2),cex.lab=1.2,cex.axis=1.2)
bar.pos <- barplot(Alpha.Group.mean,border=NA,col=col.line[Groups],ylim=c(0,y.lim),ylab="",las=1,space=0.5,xaxt="n")
mtext(side=2,text=paste(MetID,"Index",sep=" "),line=3.5,cex=1.8,adj=0.5)
mtext(side=1,text=Groups.1,at=bar.pos,line=0.5,cex=1.5)
rect(par('usr')[1],0,par('usr')[2],y.lim,border=NA)
segments(bar.pos,Alpha.Group.mean-Alpha.Group.sd,bar.pos,Alpha.Group.mean+Alpha.Group.sd)
segments(c(bar.pos-err.wid,bar.pos-err.wid),c(Alpha.Group.mean-Alpha.Group.sd,Alpha.Group.mean+Alpha.Group.sd),c(bar.pos+err.wid,bar.pos+err.wid),c(Alpha.Group.mean-Alpha.Group.sd,Alpha.Group.mean+Alpha.Group.sd),lwd=0.6)

for(m in 1:(length(Groups)-1))
{
        for(n in (m+1):length(Groups))
       {
		GG=paste(Groups[m],"-",Groups[n],sep="");
		if(Alpha.GG.Utest[[GG]]$p.value < 0.05)
		{
			par(lty=1,lwd=2,cex=2)
			lines(rep(c(bar.pos[m],bar.pos[n]-gap),each=2),c(y.lim*(1-step*(n+m))-sig.tick,y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m))-sig.tick),lwd=0.7);
			text(median(c(bar.pos[m],bar.pos[n]-gap)),y.lim*(1-step*(n+m)+0.02),labels=Alpha.GG.Asterisk[[GG]]);
		}
	}
}
dev.off()

png(paste(outDir,"/",SamID,".",MetID,".final.otu.barplot.png",sep=""),width=4*480/7,height=4*480/7, bg="transparent")
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=unique(as.character(grpInfo[,3]))
names(col.line)=unique(as.character(grpInfo[,2]))

y.lim <- max(Alpha.Group.mean+Alpha.Group.sd)*1.3;
err.wid <- 0.1
sig.tick <- 0.03
gap <- 0.008
step <- 0.03
#par(mar = c(12,14,7,6),mex=0.5,cex=1,font=1,xaxs = "i", yaxs = "i",adj=0.5,mgp=c(9,1.5,0),lty=1,lwd=2,cex.lab=2.5,cex.axis=2,tck=0.02)
par(mgp=c(4,1,0),tck=-0.02,mar=c(2,6,2,2),cex.lab=1.2,cex.axis=1.2)
bar.pos <- barplot(Alpha.Group.mean,border=NA,col=col.line[Groups],ylim=c(0,y.lim),ylab="",las=1,space=0.5,xaxt="n")
mtext(side=2,text=paste(MetID,"Index",sep=" "),line=3.5,cex=1.8,adj=0.5)
mtext(side=1,text=Groups.1,at=bar.pos,line=0.5,cex=1.5)
rect(par('usr')[1],0,par('usr')[2],y.lim,border=NA)
segments(bar.pos,Alpha.Group.mean-Alpha.Group.sd,bar.pos,Alpha.Group.mean+Alpha.Group.sd)
segments(c(bar.pos-err.wid,bar.pos-err.wid),c(Alpha.Group.mean-Alpha.Group.sd,Alpha.Group.mean+Alpha.Group.sd),c(bar.pos+err.wid,bar.pos+err.wid),c(Alpha.Group.mean-Alpha.Group.sd,Alpha.Group.mean+Alpha.Group.sd),lwd=0.6)

for(m in 1:(length(Groups)-1))
{
        for(n in (m+1):length(Groups))
       {
		GG=paste(Groups[m],"-",Groups[n],sep="");
		if(Alpha.GG.Utest[[GG]]$p.value < 0.05)
		{
			par(lty=1,lwd=2,cex=2)
			lines(rep(c(bar.pos[m],bar.pos[n]-gap),each=2),c(y.lim*(1-step*(n+m))-sig.tick,y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m))-sig.tick),lwd=0.7);
			text(median(c(bar.pos[m],bar.pos[n]-gap)),y.lim*(1-step*(n+m)+0.02),labels=Alpha.GG.Asterisk[[GG]]);
		}
	}
}
dev.off()

######### boxplot


#############boxplot

pdf(paste(outDir,"/",SamID,".",MetID,".final.otu.boxplot.pdf",sep=""),width=4,height=4)
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=unique(as.character(grpInfo[,3]))
names(col.line)=unique(as.character(grpInfo[,2]))

#par(mar = c(12,14,7,6),mex=0.5,cex=1,font=1,xaxs = "i", yaxs = "i",adj=0.5,mgp=c(9,1.5,0),lty=1,lwd=2,cex.lab=2.5,cex.axis=2,tck=0.02)
par(mgp=c(4,1,0),tck=-0.02,mar=c(2,6,2,2),cex.lab=1.2,cex.axis=1.2,bty='l')
boxplot(Alpha~groupname,border=col.line[Groups],ylab="",las=1,outline=F,xaxt="n")
mtext(side=2,text=paste(MetID,"Index",sep=" "),line=3.5,cex=1.8,adj=0.5)
xx=1:length(Groups)
mtext(side=1,text=Groups.1,at=xx,line=0.5,cex=1.5)

y.lim=par('usr')[4]

for(m in 1:(length(Groups)-1))
{
        for(n in (m+1):length(Groups))
       {
		GG=paste(Groups[m],"-",Groups[n],sep="");
		if(Alpha.GG.Utest[[GG]]$p.value < 0.05)
		{
			par(lty=1,lwd=2,cex=2)
			lines(rep(c(m,n-gap),each=2),c(y.lim*(1-step*(n+m))-sig.tick,y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m))-sig.tick),lwd=0.7);
			text(median(c(m,n-gap)),y.lim*(1-step*(n+m)+0.02),labels=Alpha.GG.Asterisk[[GG]]);
		}
	}
}
dev.off()



#############point

pdf(paste(outDir,"/",SamID,".",MetID,".final.otu.point.pdf",sep=""),width=4,height=4)
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=unique(as.character(grpInfo[,3]))
names(col.line)=unique(as.character(grpInfo[,2]))

#par(mar = c(12,14,7,6),mex=0.5,cex=1,font=1,xaxs = "i", yaxs = "i",adj=0.5,mgp=c(9,1.5,0),lty=1,lwd=2,cex.lab=2.5,cex.axis=2,tck=0.02)
par(mgp=c(4,1,0),tck=-0.02,mar=c(2,6,2,2),cex.lab=1.2,cex.axis=1.2,bty='l')
boxplot(Alpha~groupname,main="",yaxt='n',xaxt='n',medcol=col.line[Groups],boxlty = 0,ylab="",las=1,outline=F,xaxt="n",border=col.line[Groups])

u_gg=unique(groupname)
l_dd=vector('list',length=length(u_gg))
for(i in 1:length(u_gg)){
	di1=groupname %in% u_gg[i]
	l_dd[[i]]=Alpha[di1]
}
stripchart(l_dd, vertical = TRUE,method = "jitter", add = TRUE, pch = 16,col = col.line[Groups],cex=1)

mtext(side=2,text=paste(MetID,"Index",sep=" "),line=3.5,cex=1.8,adj=0.5)
xx=1:length(Groups)
mtext(side=1,text=Groups.1,at=xx,line=0.5,cex=1.5)

y.lim=par('usr')[4]

for(m in 1:(length(Groups)-1))
{
        for(n in (m+1):length(Groups))
       {
		GG=paste(Groups[m],"-",Groups[n],sep="");
		if(Alpha.GG.Utest[[GG]]$p.value < 0.05)
		{
			par(lty=1,lwd=2,cex=2)
			lines(rep(c(m,n-gap),each=2),c(y.lim*(1-step*(n+m))-sig.tick,y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m))-sig.tick),lwd=0.7);
			text(median(c(m,n-gap)),y.lim*(1-step*(n+m)+0.02),labels=Alpha.GG.Asterisk[[GG]]);
		}
	}
}
dev.off()

