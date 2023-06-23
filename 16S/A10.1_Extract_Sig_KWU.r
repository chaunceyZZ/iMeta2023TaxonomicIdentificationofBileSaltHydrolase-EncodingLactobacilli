arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}
print(arg)
#arg <-c("/home/Project/16S/HTQ_190723/Analysis/4_common_OTUs/Total/Total.otu.abundance.txt", "/home/Project/16S/HTQ_190723/bin/Total-color.txt", "/home/Project/16S/HTQ_190723/Analysis/11_Sig/Total"  ,"2" )
print(arg)
dataFile <- arg[1];
groupFile <- arg[2];
outDir <- arg[3];
grp_col <- as.numeric(arg[4]);


SamID <- strsplit(basename(dataFile),'.',fixed=T)[[1]][1];


library(vegan);
library(MASS);

grpInfo <- read.table(groupFile,header=T,check.names = FALSE, comment.char = "")
rownames(grpInfo)=as.character(grpInfo[,1])
data <- data.frame(read.table(dataFile,header=T,sep ="\t",row.names=1,check.names=F),check.names=F);
data <- data[,-ncol(data)];
Data <- data;

ii <- intersect(rownames(grpInfo),colnames(Data))
if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),ii)
	print(ii1)
}

grpInfo=grpInfo[ii,]
Data=Data[,ii]

Data_cc=colnames(Data)
groupname <- c();
for(i in 1:length(colnames(Data)))
{
	groupname <- append(groupname,as.character(grpInfo[grep(paste("^",colnames(Data)[i],"$",sep=""),grpInfo[,1]),grp_col]));
}

colnames(Data) <- groupname;
Groups <- as.character(levels(as.factor(grpInfo[,grp_col])))
Data.Group <- c();

for(i in 1:length(Groups))
{
	Data.Group[[i]] <- as.matrix(Data[,grep(paste("^",Groups[i],"$",sep=""),colnames(Data))]);
   rownames(Data.Group[[i]]) <- rownames(Data);
}

#Perform U test
if(length(Groups) >2)
{
#Extract All significant
	Data.all.Utest <- c();
	for(k in 1:length(rownames(Data))){
		Data.all.Utest[k] <- "Y";
		for(m in 1:(length(Groups)-1)){
			for(n in (m+1):length(Groups)){
				Pvalue <- wilcox.test(Data.Group[[m]][k,],Data.Group[[n]][k,])$p.value;
				if(is.na(Pvalue)){
					Data.all.Utest[k] <- "N";
				}
				else if(Pvalue > 0.05){
					Data.all.Utest[k] <- "N";
				}
			}
		}
	}
	Data.all.sig <- cbind(Data,Data.all.Utest);
	Data.all.sig <- Data.all.sig[grep("Y",Data.all.sig[,ncol(Data.all.sig)]),];
	Data.all.sig <- Data.all.sig[,-ncol(Data.all.sig)];
	colnames(Data.all.sig) <-Data_cc	
	write.table(Data.all.sig,file=paste(outDir,"/",SamID,".Utest.allSig.data.txt",sep=""),quote = F,sep="\t",col.names=NA);

#Extract at least one significant
        Data.one.Utest <- c();
        for(k in 1:length(rownames(Data))){
                Data.one.Utest[k] <- "N";
                for(m in 1:(length(Groups)-1)){
                        for(n in (m+1):length(Groups)){
                                Pvalue <- wilcox.test(Data.Group[[m]][k,],Data.Group[[n]][k,])$p.value;
                                if(is.na(Pvalue)){
                                        Data.one.Utest[k] <- "N";
                                }
                                else if(Pvalue < 0.05){
                                        Data.one.Utest[k] <- "Y";
                                }
                        }
                }
        }
        Data.one.sig <- cbind(Data,Data.one.Utest);
        Data.one.sig <- Data.one.sig[grep("Y",Data.one.sig[,ncol(Data.one.sig)]),];
        Data.one.sig <- Data.one.sig[,-ncol(Data.one.sig)];
        colnames(Data.one.sig) <- Data_cc
        write.table(Data.one.sig,file=paste(outDir,"/",SamID,".Utest.oneSig.data.txt",sep=""),quote = F,sep="\t",col.names=NA);

#Perform KW test
	Data.KWtest <- c();
	Data.KW <- as.matrix(Data);

	for(k in 1:length(rownames(Data))){
        	Data.KWtest[k] <- "N";
	        Pvalue <- kruskal.test(Data.KW[k,],factor(colnames(Data.KW)))$p.value;
	        if(is.na(Pvalue)){
	                Data.KWtest[k] <- "N";
		}
        	else if(Pvalue < 0.05){
                	Data.KWtest[k] <- "Y";
        	}
	}

	Data.kw.sig <- cbind(Data,Data.KWtest);
	Data.kw.sig <- Data.kw.sig[grep("Y",Data.kw.sig[,ncol(Data.kw.sig)]),];
	Data.kw.sig <- Data.kw.sig[,-ncol(Data.kw.sig)];
	colnames(Data.kw.sig) <- Data_cc
	write.table(Data.kw.sig,file=paste(outDir,"/",SamID,".KWtest.Sig.data.txt",sep=""),quote = F,sep="\t",col.names=NA);

}else
{
	Data.Utest <- c();
        for(k in 1:length(rownames(Data))){
                Data.Utest[k] <- "N";
                for(m in 1:(length(Groups)-1)){
                        for(n in (m+1):length(Groups)){
                                Pvalue <- wilcox.test(Data.Group[[m]][k,],Data.Group[[n]][k,])$p.value;
                                if(is.na(Pvalue)){
                                        Data.Utest[k] <- "N";
                                }
                                else if(Pvalue < 0.05){
                                        Data.Utest[k] <- "Y";
                                }
                        }
                }
        }
        Data.sig <- cbind(Data,Data.Utest);
        Data.sig <- Data.sig[grep("Y",Data.sig[,ncol(Data.sig)]),];
        Data.sig <- Data.sig[,-ncol(Data.sig)];
        colnames(Data.sig) <- Data_cc
        write.table(Data.sig,file=paste(outDir,"/",SamID,".Utest.Sig.data.txt",sep=""),quote = F,sep="\t",col.names=NA);
}


