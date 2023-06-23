arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument: Data_File Group_File Out_Dir Source_Dir Group_column1 Group_column2\n")
        quit('no')
}
print(arg)
#arg <- c("/home/wmj/02_Metagenome/16s/HTQ_190723/20190904/Analysis/1_Abundance/Total.PlyGen.abundance.txt","/home/wmj/02_Metagenome/16s/HTQ_190723/20190904/bin/Total-color.txt"  , "/home/wmj/02_Metagenome/16s/HTQ_190723/20190904/Analysis/1_Abundance" ,"/home/wmj/02_Metagenome/16s/HTQ_190723/20190904/bin")
grp_col=2
dataFile = arg[1]
groupFilie = arg[2]
outDir = arg[3]
setwd(outDir)
SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1] ####### NOTE

#data<- read.table(dataFile,header=T,row.names=1,sep='\t')
data<- read.table(dataFile,header=T,row.names=1,sep='\t',check.names = FALSE)
grpInfo <- read.table(groupFilie,header=T,check.names = FALSE, comment.char = "")
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

#Abundance <- data.frame(data[,2:ncol(data)]);
#rownames(Abundance) <- data[,1];

top <- 25
cutoff <- 0.001
dir <- dirname(dataFile)
filename <- basename(dataFile)
source_dir <- paste(arg[4],"heatmap.frq.r",sep="/")
source(source_dir)



data <- data[,as.character(grpInfo[,1])]

groupname <- c();
for(i in 1:length(colnames(data)))
{
        groupname <- append(groupname,as.character(grpInfo[grep(paste0('^',colnames(data)[i],'$'),grpInfo[,1]),2]));
}

colnames(data) <- groupname;
data_min <- min(data[data!=0]);
data[data==0]<-data_min*0.1

#data<-data*100


data4plot <- data[apply(as.matrix(data),1,function(x) any(tapply(x,grpInfo[,2],mean)>=cutoff)),]
data4plot <- cbind(data4plot,rowSums(data4plot));
data4plot <- data4plot[order(data4plot[,ncol(data4plot)],decreasing=T),];
if(nrow(data4plot) > top){
	data4plot <- data4plot[-c((top+1):nrow(data4plot)),];
}
data4plot <- data4plot[,-ncol(data4plot)];
data4plot <- data4plot[order(rownames(data4plot),decreasing=F),];
data4plot <- data4plot[,order(colnames(data4plot),decreasing=F)];

####
data4plot <- log10(data4plot)
data4plot[data4plot<(-2)]=-2;


####
Strain <- matrix(nrow=nrow(data4plot),ncol=2)
for(r in 1:nrow(Strain)){
	Strain[r,1:2] = strsplit(rownames(data4plot)[r]," ",fixed=T)[[1]]
}

Phylums <- as.character(levels(as.factor(Strain[,1])))

Strain_col <- c()
for(s in Strain[,1]){
	Strain_col <- append(Strain_col,rainbow(length(Phylums))[which(s==Phylums)])
}

#Groups <- as.character(levels(as.factor(grpInfo[,2])))
Groups <-unique(as.character(grpInfo[,grp_col]))
G_cc<-unique(as.character(grpInfo[,3]))
names(G_cc)=Groups
PopCol <- c()
for(p in colnames(data4plot)){
	id=strsplit(p,".",fixed=T)[[1]];
	ID=id[1]
	PopCol <- append(PopCol,G_cc[ID])
#	PopCol <- append(PopCol,cm.colors(length(Groups))[which(id==Groups)])
}


colors=colorRampPalette(c("navy", "white", "firebrick3"))(20)

pdf(paste(SamID,".abundance.heatmap.pdf",sep=""),height=7.5,width=12);

par(oma=c(0.1,0.1,3,0.1))
heatmap.frq(as.matrix(data4plot),ColSideColors=PopCol,Rowv=NA,Colv=NA,labRow=Strain[,2],margins=c(7,1),RowSideColors=Strain_col,verbose=F,labCol='',revC=T,cexRow=1.8,col=colors,scale='none')

par(new=T)
legend(x='bottom',fill=rainbow(length(Phylums)),Phylums,cex=1.2,horiz=T,border=NA,bty='n',xpd=T,inset=-0.1)
par(new=T,fig=c(0.14,0.3,0.86,0.97),mar=c(2.6,0.1,1,0.1),oma=c(0.1,0.1,0.1,0.1))
image(as.matrix(seq(min(data4plot),max(data4plot),length.out=length(colors))),axes=F,col=colors)
axis(1,at=seq(0,1,length.out=3),labels=round(c(seq(min(data4plot),max(data4plot),length.out=3)),0),tick=F,line=-1,cex.axis=1)
mtext(expression('log'[10]*'(relative abundance %)'),side=3,line=0,cex=1.2)


dev.off()



