arg <- commandArgs(T)
if(length(arg) < 1){
	cat("Argument: Out_Dir \n")
	quit('no')
}
print(arg)
#arg <- c("/home/Project/16S/HTQ_190723/13_betaDiv/Total/weighted_unifrac_Total.gut.sub2W2.filter_4_noChimera_rmAlignFail.otu_table_even15000.txt", "/home/Project/16S/HTQ_190723/bin/Total-color.txt"  ,"/home/Project/16S/HTQ_190723/Analysis/18_raw_beta_PCoA/Total", "2"  )


library(vegan)
library(ape) # 用于pcoa分析
library(ggplot2) #用于画图
library(ade4);
library(fpc);


outDir <- arg[3]
grp_col <- arg[4]
ff=unlist(strsplit(basename(arg[1]),'\\.'))[1]

dis <- read.table(arg[1], header=T,row.names = 1, sep = '\t')
groupInfo=read.table(arg[2],header=T,sep='\t', stringsAsFactors = FALSE, comment.char = "")
rownames(groupInfo)=as.character(groupInfo[,1])
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=unique(as.character(groupInfo[,3]))
names(col.line)=unique(as.character(groupInfo[,2]))
ii =intersect(rownames(groupInfo),rownames(dis))

if(length(ii)!=nrow(groupInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(groupInfo),ii)
	print(ii1)
}


groupInfo=groupInfo[ii,]
dis=dis[ii,ii]
pcoa =cmdscale(dis, k=2, eig=T)
Group <- unique(as.character(groupInfo[,2])) #按照目的输入样本
shape <-16 #定义点形状
color <- col.line[Group] #定义点颜色

pro1 = as.numeric(sprintf("%.3f",pcoa$eig[1]))*100
pro2 = as.numeric(sprintf("%.3f",pcoa$eig[2]))*1000
xlab=paste("PCOA1(",pro1,"%)",sep="") 
ylab=paste("PCOA2(",pro2,"%)",sep="")


x = pcoa$points
sample_names = rownames(x)
pc = as.data.frame(pcoa$points)
pc$names = sample_names
legend_title = ""
#	group = Group
#	pc$group = group
pc$group=groupInfo[rownames(pc),2]
cc=colnames(pc)
cc[1:2]=c('x','y')
colnames(pc)=cc
pdf(paste0(outDir,'/',ff,'.pcoa.pdf'))
ggplot(pc,aes(x,y)) + #用ggplot作图
	geom_point(size=3,aes(color=col.line[pc$group],shape=pc$group)) + 
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

