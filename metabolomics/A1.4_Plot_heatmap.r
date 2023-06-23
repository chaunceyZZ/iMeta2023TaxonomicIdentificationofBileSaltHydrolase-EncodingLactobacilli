arg <- commandArgs(T)
if(length(arg) != 5){
        cat("Argument: Data_File Group_File Out_Dir\n")
        quit('no')
}
print(arg)
#arg=c( "/home/wmj/02_Metagenome/Metabolome/20190904_SZW/Temp/0.0_raw_data/FHC.data.txt","/home/wmj/02_Metagenome/Metabolome/20190904_SZW/Temp/bin/FHC-grouping.info" ,"/home/wmj/02_Metagenome/Metabolome/20190904_SZW/Temp/1.5_raw_Heatmap", "2", "/home/wmj/02_Metagenome/Metabolome/20190904_SZW/Temp/bin/FHC-color.txt")
dataFile = arg[1]
groupFilie = arg[2]
outDir = arg[3]
grp_col <- as.numeric(arg[4]);

#dataFile = "FM.data.txt";
#groupFilie = "FM-grouping.info";
#outDir = "./";
#grp_col <- 2


SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1]

data<- read.table(dataFile,header=T,row.names=1,sep='\t')
grpInfo <- read.table(groupFilie,header=T)

grps <- length(grp_col)
dir <- dirname(dataFile)
filename <- basename(dataFile)
source(paste0(dirname(arg[2]),'/heatmap.frq.r'))
row.names(grpInfo)=as.character(grpInfo[,1])
data <- data[,as.character(grpInfo[,1])]
####normalize 1
#data <- scale(data, center=F,scale=T)
####normalize 2
#doit <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) - min(x, na.rm=TRUE))}
#data <- as.data.frame(lapply(data, doit))
####normalize 3
normalize <- function(x) {
    x <- sweep(x, 1, apply(x, 1, min))
    sweep(x, 1, apply(x, 1, max), "/")
}
data <- normalize(data)
data_min <- min(data[data!=0]);

data[data==0]<-data_min*0.1
#data <- log10(data);
#data[data<(-2)]=-2;

groupname <- c();
for(i in 1:length(colnames(data)))
{
        groupname <- append(groupname,as.character(grpInfo[grep(paste0('^',colnames(data)[i],"$"),as.character(grpInfo[,1])),grp_col]));
}

colnames(data) <- groupname;

####
g.c <- read.table(arg[5],header=T,row.names=1,comment.char = "")
uu=unique(g.c)
c_col=as.character(g.c[,2])
names(c_col)=as.character(g.c[,1])

PopCol <- c_col[colnames(data)]

colors=colorRampPalette(c("green", "black", "firebrick3"))(256)

Hclust <- c("complete","average","median","centroid");
Corrlations <- c("pearson", "kendall", "spearman");
for (CID in Corrlations){
	for(HID in Hclust){
		pdf(paste(outDir,"/",SamID,".",CID,".",HID,".heatmap.pdf",sep=""),height=7.5,width=12);

		distCor <- function(x) as.dist(1-cor(t(x),method=CID));
		hclustMet <-function(x) hclust(x,method=HID);

		par(oma=c(0.1,0.1,3,0.1))
		heatmap.frq(as.matrix(data),ColSideColors=PopCol,distfun = distCor,hclust = hclustMet, labRow=NA,margins=c(7,1),verbose=F,labCol=groupname,revC=T,cexRow=1.8,col=colors,scale='none');

		par(new=T)
		par(new=T,fig=c(0.14,0.3,0.86,0.97),mar=c(2.1,2.1,1,0.1),oma=c(0.1,0.1,0.1,0.1))
		image(as.matrix(seq(min(data),max(data),length.out=length(colors))),axes=F,col=colors)
		axis(1,at=seq(0,1,length.out=3),labels=round(c(seq(min(data),max(data),length.out=3)),0),tick=F,line=-1,cex.axis=1)
        #mtext(expression('log'[10]*'(Abundance)'),side=3,line=0,cex=1.2)
		mtext('Normalized abundance',side=3,line=0,cex=1.2)


		##	legend
		par(new=T,fig=c(0.05,0.3,0.74,0.85),mar=c(2.1,2.1,1,0.1),oma=c(0.1,0.1,0.1,0.1))
		plot(0,0,col="white",xlab="",ylab="",axes=F)
		legend("left",fill=unique(PopCol),legend=unique(g.c[,1]),horiz=T,xpd=T,bty="n")
		dev.off()
	}
}



