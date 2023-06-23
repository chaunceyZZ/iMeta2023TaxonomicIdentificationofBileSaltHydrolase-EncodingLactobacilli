arg <- commandArgs(T)
if(length(arg) != 4){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}
print(arg)
#arg <-c("Total.cog.L3.abundance.txt","tt.4lefse.txt","Total-grouping.info",2)
data <- read.table(arg[1],sep="\t",header=T)
grp <- read.table(arg[3],sep="\t",header=T, comment.char = "")
rownames(grp)<-as.character(grp[,1])
ii <-intersect(rownames(grp),colnames(data))
if(length(ii)!=nrow(grp)){
	print('Data Lost')
	ii1=setdiff(rownames(grp),ii)
	print(ii1)
}

gr <-grp[ii,]
dd <-data[,c("ID",ii)]
cc <- as.character(c("Class",as.character(gr[,2])))
dd1 <-rbind(cc,colnames(dd),as.matrix(dd))
write.table(dd1,arg[2],sep="\t",quote=F,row.names=F,col.names=F)
