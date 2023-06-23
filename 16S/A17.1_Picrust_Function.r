arg <- commandArgs(T)
if(length(arg) < 1){
	cat("Argument: Out_Dir \n")
	quit('no')
}
print(arg)
#arg <- c("/home/wmj/02_Metagenome/16s/JZM_20190119/20190325/Analysis/20_PICRUSt_Function/Total.catagrized.kegg.L2_L2.txt", "/home/wmj/02_Metagenome/16s/JZM_20190119/20190325/Analysis/20_PICRUSt_Function/Total.kegg.L3.res" ,"/home/wmj/02_Metagenome/16s/JZM_20190119/20190325/bin/Total-grouping.info" , "/home/wmj/02_Metagenome/16s/JZM_20190119/20190325/Analysis/20_PICRUSt_Function"    )

gene=read.table(arg[1],sep='\t',header=T,row.names=1)
kegg=read.table(arg[2],row.names=1,sep="\t")
group=read.table(arg[3],header=T, comment.char = "")
rownames(group)=as.character(group[,1])
#col.line=c("#DC143C","#4169E1","#2E8B57","#9932CC","#FF8C00","#FFFF00","#8B4513","#FF69B4","#808080")
col.line=unique(as.character(group[,3]))
names(col.line)=unique(as.character(group[,2]))
################	gene
rr=rownames(gene)
r1=c()
for(i in 1:length(rr)){
	r1=c(r1,unlist(strsplit(as.character(rr[i]),';'))[2])
}
r2=gsub(' ','',r1)
###### kegg
k1=is.na(kegg[,3])
kk=kegg[!k1,]
#######
ii=intersect(r2,rownames(kk))
rownames(gene)=r2
gene1=gene[ii,]
gr=unique(as.character(group[,2]))
names(r1)=r2
for(i in 1:length(ii)){
	d1=gene1[i,]
	ss=sum(d1)
	d2=d1/ss
	tt=order(names(d2))
	t1=names(d2)[tt]
	d3=d2[t1]
	gg=as.character(group[names(d3),2])

	rr=matrix(as.numeric(d3),byrow=F,ncol=length(gr))
	colnames(rr)=gr
	ff=paste0(arg[4],'/',ii[i],'.pdf')
	pdf(ff)
	par(bty='l')
#	xx=boxplot(as.numeric(d3)~gg,ylab='Proportion of sequences (%)',outline=F,main=r1[ii[i]],col=col.line[1:length(gr)],border=col.line[1:length(gr)])
	xx=boxplot(rr,ylab='Proportion of sequences (%)',outline=F,main=r1[ii[i]],col=col.line[gr],border=col.line[gr])

	gap=par('usr')[4]
	for(m in 1:(length(gr)-1)){
		for(n in (m+1):length(gr)){
			g1=grep(paste0("^",gr[m],'$'),as.character(group[,2]))
			g2=grep(paste0('^',gr[n],'$'),as.character(group[,2]))
			
			t1=d1[,rownames(group)[g1]]
			t2=d1[,rownames(group)[g2]]

			p.t=t.test(t1,t2,paired = FALSE)$p.value
			if(p.t<0.05){ f1<- '*'}
			if(p.t<0.01){f1 <- '**'}
			if(p.t<0.001){f1 <- '***'}
			if(p.t <0.05){
				nn=0.001*m+0.001*n
				lines(x=c(m+0.02,n-0.02),y=c(gap+nn,gap+nn),lwd=0.7,xpd=T);
				text(x=median(c(m+0.02,n-0.02)),y=gap+nn+0.001,labels=f1,xpd=T,adj=0.5);
			}
		}
	}




	dev.off()
}


