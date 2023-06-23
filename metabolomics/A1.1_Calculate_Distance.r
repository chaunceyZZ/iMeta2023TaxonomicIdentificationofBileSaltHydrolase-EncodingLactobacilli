arg <- commandArgs(T)
if(length(arg) != 2){
        cat("Argument: Data_File Out_Dir\n")
        quit('no')
}
print(arg)
#arg=c("/home/wmj/02_Metagenome/Metabolome/HFQ/外包/20181224_HFQ/0.0_raw_data/NC.data.txt","/home/wmj/02_Metagenome/Metabolome/HFQ/外包/20181224_HFQ/1.2_raw_Distance" )
dataFile <- arg[1];
outDir <- arg[2];

#dataFile <- "total.Utest.oneSig.data.txt";
#outDir <- "./";

library(ecodist);
FileName <- strsplit(basename(dataFile),'.',fixed=T)[[1]];
SamID <- FileName[1];

#Data <- t(read.table(dataFile,header=T,row.names = 1));
dd=read.table(dataFile,header=T,row.names = 1)


Corrlations <- c("pearson", "kendall", "spearman");
for (CID in Corrlations){
#Type <-1- cor(t(Data),method=CID);
Type <-1- cor(dd,method=CID);
write.table(Type,file=paste(outDir,"/",SamID,".",CID,".txt",sep=""),quote = F,sep="\t",col.names=NA);
}

