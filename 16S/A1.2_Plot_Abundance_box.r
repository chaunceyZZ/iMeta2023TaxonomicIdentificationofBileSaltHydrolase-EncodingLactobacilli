arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument: Data_File Group_File Out_Dir\n")
        quit('no')
}
print(arg)
#arg <- c("/home/wmj/02_Metagenome/16s/HTQ_190723/Analysis/1_Abundance/Total.Genus.abundance.txt", "/home/wmj/02_Metagenome/16s/HTQ_190723/bin/Total-color.txt" ,"/home/wmj/02_Metagenome/16s/HTQ_190723/Analysis/1_Abundance"    ,"15" )


dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
setwd(outDir)

SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1] ####### NOTE


top <- as.numeric(arg[4])

data <- read.table(dataFile,header=T,row.names=1,check.names = FALSE);
grpInfo <- read.table(groupFile,header=T,check.names = FALSE, comment.char = "")
rownames(grpInfo) <- grpInfo[,1]
grp <- basename(outDir)
#nn <- grep(paste("^",grp,".*",sep=""),colnames(grpInfo))
#nn1 <- colnames(grpInfo)[nn]
#nn2 <- strsplit(nn1,".",fixed=T)[[1]][2]
#nn3 <- grep(paste("^",nn2,".*",sep=""),colnames(grpInfo))
ii <- intersect(rownames(grpInfo),colnames(data))
if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),colnames(data))
	print(ii1)
}
grpInfo <- grpInfo[ii,]
data <-data[,ii]


grp_col=2
Abundance <-data;

groupname <- c();
for(i in 1:length(colnames(Abundance)))
{
        groupname <- append(groupname,as.character(grpInfo[grep(paste0("^",colnames(Abundance)[i],'$'),grpInfo[,1]),grp_col]));
}

colnames(Abundance) <- groupname;
#Groups <- as.character(levels(as.factor(grpInfo[,grp_col])))
Groups <-unique(as.character(grpInfo[,grp_col]))
G_cc<-unique(as.character(grpInfo[,3]))
names(G_cc)=Groups
Abundance.Group <- c();
for(i in 1:length(Groups))
{
        Abundance.Group[[i]] <- as.matrix(Abundance[,grep(paste("^",Groups[i],"$",sep=""),colnames(Abundance))]);
	rownames(Abundance.Group[[i]]) <- rownames(Abundance);
}

#######Calculate the abundance of each type and boxplot the distribution
pdf(paste(SamID,".abundance.boxplot.pdf",sep=""),width=14,height=7);

#par(mfrow=c(1,1),mar=c(30,3,3,3));

Type.abundance=as.matrix(data.frame(Abundance[,1:ncol(Abundance)],sum=rowSums(Abundance[,1:ncol(Abundance)])));
Type.abundance=Type.abundance[order(Type.abundance[,ncol(Type.abundance)],decreasing=T),];	# sorted by the total abundance
Type.abundance=Type.abundance[,-ncol(Type.abundance)];
colnames(Type.abundance) <- colnames(Abundance)
Type.abundance.Group <- c()
for(i in 1:length(Groups))
{
        Type.abundance.Group[[i]] <- as.matrix(Type.abundance[,grep(paste("^",Groups[i],"$",sep=""),colnames(Type.abundance))]);
        rownames(Type.abundance.Group[[i]]) <- rownames(Type.abundance);
}
Top.abundance <-  Type.abundance;
if(nrow(Top.abundance)>top){															
# get top types
	Top.abundance=Top.abundance[-c((top+1):nrow(Top.abundance)),];
}

Top.abundance.Group  <- c();
for(i in 1:length(Groups))
{
        Top.abundance.Group[[i]] <- as.matrix(Top.abundance[,grep(paste("^",Groups[i],"$",sep=""),colnames(Top.abundance))]);
        rownames(Top.abundance.Group[[i]]) <- rownames(Top.abundance);
}


Top.group=matrix(nrow=nrow(Top.abundance),ncol=ncol(Top.abundance));
Top.color=matrix(nrow=nrow(Top.abundance),ncol=ncol(Top.abundance));
Top.name=rownames(Top.abundance);

for(m in 1:nrow(Top.abundance))
{
	for(n in 1:ncol(Top.abundance))
	{
		Top.group[m,n]= (m-1)*length(Groups)+ which(colnames(Top.abundance)[n]==Groups);
	}
}


####perform wilcox U test to compare the abundance of each type among three populations####
Results=matrix(nrow=nrow(Type.abundance),ncol=3);
rownames(Results)=rownames(Type.abundance);
colnames(Results)=c("Aov.p","Aov.pBH","Aov.pbonferroni");
Top.significant <- rownames(Top.abundance);
Results[,1] <- apply(Type.abundance,1,function(x) summary(aov(x~groupname))[[1]][,5][1])
Results[,2] <- p.adjust(Results[,1],'BH')
Results[,3] <- p.adjust(Results[,1],'bonferroni')

top_top_n = 0;
for(r in 1:nrow(Top.abundance)){
if(mean(Top.abundance[r,])>=0.05){top_top_n = r}
else{break}
}

Top.abundance.perc = Top.abundance*100;
num_groups <- length(Groups);
rect_bot <- seq(0.5,max(Top.group),2*num_groups)
rect_top <- seq(num_groups+0.5,max(Top.group)+0.5,2*num_groups)


par(mfrow=c(1,1),mar=c(16.3,4,4,4));
boxplot(Top.abundance.perc[1:top_top_n,]~Top.group[1:top_top_n,],col=G_cc[Groups],xaxt="n",xlim=c(1,max(Top.group)),outline=F,ylab='',cex.axis=1.3,xlab='');
text(x=seq((num_groups-1),(num_groups*top-1),num_groups),y=par("usr")[3]-0.5,labels=Top.name,xpd=T,srt=45,pos=2,font=3,cex=1.1,offset=0)

rect(rect_bot[1:ceiling(top_top_n/2)],-10,rect_top[1:ceiling(top_top_n/2)],100,col='grey90',lty=0)
boxplot(Top.abundance.perc[1:top_top_n,]~Top.group[1:top_top_n,],col=G_cc[Groups],xaxt="n",add=T,outline=F,yaxt='n',xlab='');

par(new=T);

boxplot(Top.abundance.perc[(top_top_n+1):nrow(Top.abundance.perc),]~Top.group[(top_top_n+1):nrow(Top.abundance.perc),],col=G_cc[Groups],xaxt="n",xlim=c(1,max(Top.group)),at=min(Top.group[(top_top_n+1):nrow(Top.abundance.perc),]):max(Top.group[(top_top_n+1):nrow(Top.abundance.perc),]),yaxt='n',outline=F,ylab='',xlab='');
rect(rect_bot[(ceiling(top_top_n/2)+1):length(rect_bot)],-10,rect_top[(ceiling(top_top_n/2)+1):length(rect_top )],100,col='grey90',lty=0)
boxplot(Top.abundance.perc[(top_top_n+1):nrow(Top.abundance.perc),]~Top.group[(top_top_n+1):nrow(Top.abundance.perc),],col=G_cc[Groups],xlim=c(1,max(Top.group)),at=min(Top.group[(top_top_n+1):nrow(Top.abundance.perc),]):max(Top.group[(top_top_n+1):nrow(Top.abundance.perc),]),add=T,xaxt='n',yaxt='n',outline=F,xlab='');
abline(v=max(Top.group[1:top_top_n,])+0.5,lty=2,lwd=3,col="red")
legend(x='topright',fill=G_cc[Groups],bty='n',legend=Groups,cex=1.2)
axis(4,cex.axis=1.3)
mtext(side=2,line=2,'relative abundance (%)',cex=1.5)
mtext(side=4,line=2,'relative abundance (%)',cex=1.5)
vv <-max(Top.group[1:top_top_n,])+0.5
mtext(side=3,line=2,at=c(vv-0.5,vv+0.5),text=c(expression('R.A.'>='5%'),expression('R.A.'<'5%')),cex=1.2,adj=c(1,0))
arrows(x0=vv-0.5,y0=par("usr")[4]+1,x1=vv-4.5,y1=par("usr")[4]+1,col="black",lwd=2,length=0.05,code=2,angle=45,xpd=T)
arrows(x0=vv+0.5,y0=par("usr")[4]+1,x1=vv+4.5,y1=par("usr")[4]+1,col="black",lwd=2,length=0.05,code=2,angle=45,xpd=T)


#mtext(side=3,line=1,at=c(mean(rect_bot[1:ceiling(top_top_n/2)]),mean(rect_bot[(ceiling(top_top_n/2)+1):length(rect_bot)])),text=c(expression('R.A.'>='5%'),expression('R.A.'<'5%')),cex=1.2)

for(i in 1:nrow(Top.abundance)){
	Top.significant[i] <- '   ';
	if(Results[i,'Aov.p']<0.05){Top.significant[i] <- ' * '}
	if(Results[i,'Aov.p']<0.01){Top.significant[i] <- ' ** '}
	if(Results[i,'Aov.p']<0.001){Top.significant[i] <- '***'}
}

axis(3,labels=Top.significant,at=seq((num_groups-1),(num_groups*top-1),num_groups),tick = F,line = -1.2);
dev.off()

################
png(paste(SamID,".abundance.boxplot.png",sep=""),width=14*480/7,height=7*480/7, bg="transparent")
par(mfrow=c(1,1),mar=c(16.3,4,4,4));
boxplot(Top.abundance.perc[1:top_top_n,]~Top.group[1:top_top_n,],col=G_cc[Groups],xaxt="n",xlim=c(1,max(Top.group)),outline=F,ylab='',cex.axis=1.3,xlab='');
text(x=seq((num_groups-1),(num_groups*top-1),num_groups),y=par("usr")[3]-0.5,labels=Top.name,xpd=T,srt=45,pos=2,font=3,cex=1.1,offset=0)

rect(rect_bot[1:ceiling(top_top_n/2)],-10,rect_top[1:ceiling(top_top_n/2)],100,col='grey90',lty=0)
boxplot(Top.abundance.perc[1:top_top_n,]~Top.group[1:top_top_n,],col=G_cc[Groups],xaxt="n",add=T,outline=F,yaxt='n',xlab='');

par(new=T);

boxplot(Top.abundance.perc[(top_top_n+1):nrow(Top.abundance.perc),]~Top.group[(top_top_n+1):nrow(Top.abundance.perc),],col=G_cc[Groups],xaxt="n",xlim=c(1,max(Top.group)),at=min(Top.group[(top_top_n+1):nrow(Top.abundance.perc),]):max(Top.group[(top_top_n+1):nrow(Top.abundance.perc),]),yaxt='n',outline=F,ylab='',xlab='');
rect(rect_bot[(ceiling(top_top_n/2)+1):length(rect_bot)],-10,rect_top[(ceiling(top_top_n/2)+1):length(rect_top )],100,col='grey90',lty=0)
boxplot(Top.abundance.perc[(top_top_n+1):nrow(Top.abundance.perc),]~Top.group[(top_top_n+1):nrow(Top.abundance.perc),],col=G_cc[Groups],xlim=c(1,max(Top.group)),at=min(Top.group[(top_top_n+1):nrow(Top.abundance.perc),]):max(Top.group[(top_top_n+1):nrow(Top.abundance.perc),]),add=T,xaxt='n',yaxt='n',outline=F,xlab='');
abline(v=max(Top.group[1:top_top_n,])+0.5,lty=2,lwd=3,col="red")
legend(x='topright',fill=G_cc[Groups],bty='n',legend=Groups,cex=1.2)
axis(4,cex.axis=1.3)
mtext(side=2,line=2,'relative abundance (%)',cex=1.5)
mtext(side=4,line=2,'relative abundance (%)',cex=1.5)
vv <-max(Top.group[1:top_top_n,])+0.5
mtext(side=3,line=2,at=c(vv-0.5,vv+0.5),text=c(expression('R.A.'>='5%'),expression('R.A.'<'5%')),cex=1.2,adj=c(1,0))
arrows(x0=vv-0.5,y0=par("usr")[4]+1,x1=vv-4.5,y1=par("usr")[4]+1,col="black",lwd=2,length=0.05,code=2,angle=45,xpd=T)
arrows(x0=vv+0.5,y0=par("usr")[4]+1,x1=vv+4.5,y1=par("usr")[4]+1,col="black",lwd=2,length=0.05,code=2,angle=45,xpd=T)


#mtext(side=3,line=1,at=c(mean(rect_bot[1:ceiling(top_top_n/2)]),mean(rect_bot[(ceiling(top_top_n/2)+1):length(rect_bot)])),text=c(expression('R.A.'>='5%'),expression('R.A.'<'5%')),cex=1.2)

for(i in 1:nrow(Top.abundance)){
	Top.significant[i] <- '   ';
	if(Results[i,'Aov.p']<0.05){Top.significant[i] <- ' * '}
	if(Results[i,'Aov.p']<0.01){Top.significant[i] <- ' ** '}
	if(Results[i,'Aov.p']<0.001){Top.significant[i] <- '***'}
}

axis(3,labels=Top.significant,at=seq((num_groups-1),(num_groups*top-1),num_groups),tick = F,line = -1.2);

dev.off()

