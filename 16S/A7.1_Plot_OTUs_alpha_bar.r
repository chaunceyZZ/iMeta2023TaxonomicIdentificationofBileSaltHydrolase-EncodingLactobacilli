arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}
print(arg)
#arg <- c( "/home/wmj/02_Metagenome/16s/HTQ_190723/14_Statistics/alpha/Total/Total.final.otu.chao1.txt","/home/wmj/02_Metagenome/16s/HTQ_190723/bin/Total-color.txt", "/home/wmj/02_Metagenome/16s/HTQ_190723/Analysis/7_OTUs_alpha/Total"   , "2"  )

dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);

grpInfo <- read.table(groupFile,header=T,check.names = FALSE, comment.char = "")
rownames(grpInfo)<-as.character(grpInfo[,1])
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
	groupname <- append(groupname,as.character(grpInfo[grep(paste("^",rownames(Alpha)[i],"$",sep=""),as.character(grpInfo[,1])),grp_col]));
}

rownames(Alpha) <- groupname;
#Groups <- as.character(levels(as.factor(grpInfo[,grp_col])))
Groups <- as.character(unique(grpInfo[,grp_col]))
G_cc=as.character(unique(grpInfo[,3]))
names(G_cc)=Groups
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
col.line=G_cc
y.lim <- max(Alpha.Group.mean+Alpha.Group.sd)*1.3;
err.wid <- 0.1
sig.tick <- 0.03
gap <- 0.008
step <- 0.03

bar.pos <- barplot(Alpha.Group.mean,border=NA,col=col.line[Groups],names.arg="",ylim=c(0,y.lim),ylab=paste(MetID,"Index",sep=" "),cex.names=1.0,cex.lab=1.0,cex=1.2)
#text(x=seq((num_groups-1),(num_groups*top-1),num_groups),y=par("usr")[3]-0.5,labels=Top.name,xpd=T,srt=45,pos=2,font=3,cex=1.1,offset=0)

text(x=bar.pos,y=par("usr")[3]-0.5,labels=Groups,xpd=T,srt=45,pos=2,font=3,cex=1,offset=0)
#mtext(at=bar.pos,text=Groups,srt=-45)
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
			lines(rep(c(bar.pos[m],bar.pos[n]-gap),each=2),c(y.lim*(1-step*(n+m))-sig.tick,y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m)),y.lim*(1-step*(n+m))-sig.tick),lwd=0.7);
			text(median(c(bar.pos[m],bar.pos[n]-gap)),y.lim*(1-step*(n+m)+0.02),labels=Alpha.GG.Asterisk[[GG]]);
		}
	}
}
		

dev.off()
