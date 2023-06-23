File <- Sys.glob("total*.dat")     #文件名称
for (f in 1:length(File))
{
	file = File[f]
	data<-read.table(file,header=T)
	count<-data.frame(table(data[,1]))
	data.new<-as.matrix(data)
    name<-paste(count[grep(2,count[,2]),1],".2",sep="")
    data.new[grep(count[grep(2,count[,2]),1],data.new[,1])[2],1]<-name
	write.table(data.new,paste("change.",file,sep=""),quote=F,sep="\t",row.names=F)
}
#data[grep(count[grep(2,count[,2]),1],data[,1])[2],1]<-paste(data[grep(count[grep(2,count[,2]),1],data[,1])[2],1],".2",sep="")
