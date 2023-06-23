arg <- commandArgs(T)
if(length(arg) < 1){
	cat("Argument: Out_Dir \n")
	quit('no')
}
print(arg)
#arg <- c("/home/wmj/02_Metagenome/16s/JZM_20190119/20190325/Analysis/12_Sig_Distance/Total/Total.KWtest.Sig.difference.txt","/home/wmj/02_Metagenome/16s/JZM_20190119/20190325/bin/Total-grouping.info"   ,"/home/wmj/02_Metagenome/16s/JZM_20190119/20190325/Analysis/19_Sig_PCoA/Total"     )


library(vegan)
library(ape) # 用于pcoa分析
library(ggplot2) #用于画图
library(ade4);
library(fpc);
outDir <- arg[3];
ff=unlist(strsplit(basename(arg[1]),'\\.txt'))[1]
dis <- read.table(arg[1], header=T,row.names = 1, sep = '\t')
group=read.table(arg[2],header=T,sep='\t', stringsAsFactors = FALSE, comment.char = "")
rownames(group)=as.character(group[,1])
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=unique(as.character(group[,3]))
names(col.line)=unique(as.character(group[,2]))
ii =intersect(rownames(group),rownames(dis))
if(length(ii)!=nrow(group)){
	print('Data Lost')
	ii1=setdiff(rownames(group),ii)
	print(ii1)
}



group=group[ii,]
dis=dis[ii,ii]

Group <- unique(as.character(group[,2])) #按照目的输入样本
shape <-16 #定义点形状
color <- col.line[Group] #定义点颜色
PCOA <- pcoa(dis, correction="none", rn=NULL) #利用PCOA()指令做pcoa分析
result <-PCOA$values[,"Relative_eig"]
pro1 = as.numeric(sprintf("%.3f",result[1]))*100
pro2 = as.numeric(sprintf("%.3f",result[2]))*100
x = PCOA$vectors
sample_names = rownames(x)
pc = as.data.frame(PCOA$vectors)
pc$names = sample_names
legend_title = ""
#group = Group
#pc$group = group
pc$group=group[rownames(pc),2]
xlab=paste("PCOA1(",pro1,"%)",sep="") 
ylab=paste("PCOA2(",pro2,"%)",sep="")
pdf(paste0(outDir,'/',ff,'.pcoa.pdf'))
ggplot(pc,aes(Axis.1,Axis.2)) + #用ggplot作图
  geom_point(size=3,aes(color=group,shape=group)) + 
#  geom_text(aes(label=names),size=4,vjust=-1) +
  labs(x=xlab,y=ylab,title="PCOA",color=legend_title,shape=legend_title) + 
  geom_hline(yintercept=0,linetype=4,color="grey") + 
  geom_vline(xintercept=0,linetype=4,color="grey") + 
  theme_bw()
dev.off()


pdf(paste0(outDir,'/',ff,'.pcoa.circle.pdf'))

plot(pc[,1],pc[,2],col = color,pch=shape,xlab=xlab,ylab=ylab,cex=0.9);

s.class(pc[,1:2], fac=as.factor(pc$group),grid=F, addaxes=F,axesell =T,label=Group,col=col.line[Group],pch=shape,add.plot=T);

dev.off()

pdf(paste0(outDir,'/',ff,'.pcoa.no.circle.pdf'))

plot(pc[,1],pc[,2],col = color,pch=shape,xlab=xlab,ylab=ylab,cex=0.9);

legend(x='topright',fill=color,bty='n',legend=Group,cex=1.2)
dev.off()

