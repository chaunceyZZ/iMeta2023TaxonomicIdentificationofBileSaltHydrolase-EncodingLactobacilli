library("muma")


arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}
print(arg)
#arg <-c(  "/home/Project/Metabolome/HFQ/HFQ_LOESS/HILIC/0.0_raw_data/FABAM.data.txt", "/home/Project/Metabolome/HFQ/HFQ_LOESS/HILIC/bin/FABAM-grouping.info","/home/Project/Metabolome/HFQ/HFQ_LOESS/HILIC/3.2_OPLS_DA" ,"2" ,"/home/Project/Metabolome/HFQ/HFQ_LOESS/HILIC" )
dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);

#dataFile <- "total.data.txt";
#groupFile <- "total-grouping.info";
#outDir <- "./";
#grp_col <- 2;


grpInfo <- read.table(groupFile,header=T);
FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- FileName[1];
groupname <- c();
Data <- read.table(dataFile,header=T,sep="\t");
nn <- as.character(Data[,1])
dd <- as.matrix(Data[,-1])
rownames(dd) <-nn
Data <-dd
Data.tm <- t(as.matrix(Data)); 
for(i in 1:length(colnames(Data)))
{
	groupname <- append(groupname,as.numeric(grpInfo[grep(paste("^",colnames(Data)[i],"$",sep=""),grpInfo[,1]),grp_col]));
}
Data.muma <- cbind(groupname,Data.tm);
outFile <- paste(SamID,".muma.csv",sep="");
write.csv(Data.muma,file=outFile);

#source("/home/Metagenome-2/Metabolome/ZX-1/bin/explore.data.lj.R")
#source("/home/Metagenome-2/Metabolome/ZX-1/bin/oplsda.lj.R")
source(paste(arg[5],"/bin/explore.data.lj.R",sep=""))
source(paste(arg[5],"/bin/oplsda.lj.R",sep=""))
explore.data(outFile,scaling="P", scal = TRUE, normalize = TRUE, imputation = FALSE, imput=c("NA"));

oplsda("P");


