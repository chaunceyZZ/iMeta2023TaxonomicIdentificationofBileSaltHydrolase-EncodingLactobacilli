arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}
print(arg)
#arg <- c( "/home/wmj/02_Metagenome/Metabolome/20190904_SZW/20190905/0.0_raw_data/P.data.txt","/home/wmj/02_Metagenome/Metabolome/20190904_SZW/20190905/bin/P-grouping.info" , "/home/wmj/02_Metagenome/Metabolome/20190904_SZW/20190905/2.1_Sig_KWU" , "2" )


dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);

#dataFile <- "../0.0_raw_data/CY.data.txt";
#groupFile <- "CY-grouping.info";
#outDir <- "../2.1_Sig_KWU";
#grp_col <- 2;

SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1];

library(vegan);
library(MASS);

grpInfo <- read.table(groupFile,header=T);
data <-read.table(dataFile,header=T,sep ="\t");
nn <- as.character(data[,1])
dd <- as.matrix(data[,-1])
rownames(dd) <-nn
Data <- dd;

groupname <- c();
for(i in 1:length(colnames(Data)))
{
	groupname <- append(groupname,as.character(grpInfo[grep(paste("^",colnames(Data)[i],"$",sep=""),as.character(grpInfo[,1])),grp_col]));
}

colnames(Data) <- groupname;

Groups <- unique(groupname)
Data.Group <- c();
for(i in 1:length(Groups))
{
	Data.Group[[i]] <- as.matrix(Data[,grep(paste("^",Groups[i],"$",sep=""),colnames(Data))]);
        rownames(Data.Group[[i]]) <- rownames(Data);
}

#Perform U test
if(length(Groups) >2)
{

#Perform KW test
	Data.KWtest <- c();
	Data.KW <- as.matrix(Data);
	
	for(k in 1:length(rownames(Data))){
        #Data.KWtest[k] <- "N";
        #Pvalue <- kruskal.test(Data.KW[k,],factor(colnames(Data.KW)))$p.value;
        #if(is.na(Pvalue)){
        #Data.KWtest[k] <- "N";
        #}
        #else if(Pvalue < 0.05){
        #Data.KWtest[k] <- "Y";
        #}

        	Data.KWtest[k] <- kruskal.test(Data.KW[k,],factor(colnames(Data.KW),levels=unique(colnames(Data.KW)),labels=unique(colnames(Data.KW))))$p.value
	}
	Data.KWtest[is.na(Data.KWtest)] <- 1
	names(Data.KWtest) <- rownames(Data)
	write.table(Data.KWtest,file=paste(outDir,"/",SamID,".All.p.value.txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)

	Data.KWtest <- p.adjust(Data.KWtest,method='fdr')
	names(Data.KWtest) <- rownames(Data)
	write.table(Data.KWtest,file=paste(outDir,"/",SamID,".All.p.adjust.txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)
	t.kw <- Data.KWtest[which(Data.KWtest < 0.05)]
	write.table(t.kw,file=paste(outDir,"/",SamID,".Sig.p.adjust.txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)
#Data.kw.sig <- cbind(Data,Data.KWtest);
#	Data.kw.sig <- Data.kw.sig[grep("Y",Data.kw.sig[,ncol(Data.kw.sig)]),];
#	Data.kw.sig <- Data.kw.sig[,-ncol(Data.kw.sig)];
    Data.kw.sig <- Data[Data.KWtest < 0.05,]
	colnames(Data.kw.sig) <- colnames(data)[-1]
    write.table(Data.kw.sig,file=paste(outDir,"/",SamID,".KWtest.Sig.data.txt",sep=""),quote = F,sep="\t",col.names=NA);

}else
{
	Data.Utest <- c();
        for(k in 1:length(rownames(Data))){
            #Data.Utest[k] <- "N";
                for(m in 1:(length(Groups)-1)){
                        for(n in (m+1):length(Groups)){
                            Pvalue <- wilcox.test(Data.Group[[m]][k,],Data.Group[[n]][k,])$p.value;
                            Data.Utest <- c(Data.Utest, Pvalue)
                            #   if(is.na(Pvalue)){
                            #           Data.Utest[k] <- "N";
                            #   }
                            #   else if(Pvalue < 0.05){
                            #           Data.Utest[k] <- "Y";
                            #   }
                        }
                }
        }
        Data.Utest[is.na(Data.Utest)] <- 1
	names(Data.Utest) <- rownames(Data)
	write.table(Data.Utest,file=paste(outDir,"/",SamID,".All.p.value.txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)

        Data.Utest <- p.adjust(Data.Utest,method='fdr')
	names(Data.Utest) <- rownames(Data)
	write.table(Data.Utest,file=paste(outDir,"/",SamID,".All.p.adjust.txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)
	t.kw <- Data.Utest[which(Data.Utest < 0.05)]
	write.table(t.kw,file=paste(outDir,"/",SamID,".Sig.p.adjust.txt",sep=""),sep="\t",quote=F,row.names=T,col.names=NA)
        #Data.sig <- cbind(Data,Data.Utest);
        #Data.sig <- Data.sig[grep("Y",Data.sig[,ncol(Data.sig)]),];
        #Data.sig <- Data.sig[,-ncol(Data.sig)];
        #colnames(Data.sig) <- colnames(data)
        Data.sig <- Data[Data.Utest < 0.05,]
	if(class(Data.sig) == "numeric"){
		Data.sig <- t(as.matrix(Data.sig))
		rownames(Data.sig) <- rownames(Data)[Data.Utest < 0.05]
		colnames(Data.sig) <- colnames(data)[-1]
	} else {
		colnames(Data.sig) <- colnames(data)[-1]
	}	
        write.table(Data.sig,file=paste(outDir,"/",SamID,".Utest.Sig.data.txt",sep=""),quote = F,sep="\t",col.names=NA);
}


