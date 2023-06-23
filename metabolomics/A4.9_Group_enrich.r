arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}

dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);

#dataFile <- "total.data.txt";
#groupFile <- "total-grouping.info";
#outDir <- "./";
#grp_col <- 2;

SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1];
MID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][2];
LID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][3];

top <- 20

data <- read.table(dataFile,header=T);
grpInfo <- read.table(groupFile,header=T);
Abundance <- data.frame(data[,2:ncol(data)]);
rownames(Abundance) <- data[,1];

groupname <- c();
for(i in 1:length(colnames(Abundance)))
{
        groupname <- append(groupname,as.character(grpInfo[grep(paste0('^',colnames(Abundance)[i],'$'),as.character(grpInfo[,1])),grp_col]));
}

colnames(Abundance) <- groupname;
Groups <- unique(as.character(grpInfo[,grp_col]))
Abundance.Group <- c();
Average.Group <- c();
for(i in 1:length(Groups))
{
        Abundance.Group[[i]] <- as.matrix(Abundance[,grep(paste0('^',Groups[i],'$'),colnames(Abundance))]);
	Average.Group <- cbind(Average.Group,apply(Abundance.Group[[i]],1,mean));
	rownames(Abundance.Group[[i]]) <- rownames(Abundance);
}

colnames(Average.Group) <- Groups;

maxname <- function(x) names(which.max(x));
Max.Group <- apply(Average.Group,1,maxname);
write.table(Max.Group,file=paste(outDir,"/",SamID,".",MID,".",LID,".VIP.enriched.metabolite.dat",sep=""),quote = F,sep="\t",col.names=NA);
	
