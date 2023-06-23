
arg <- commandArgs(T)
if(length(arg) != 3){
	cat("Argument: dataFile GroupInfo Out_Dir Group_column Y C\n")
	quit('no')
}

print(arg)
#arg <-c( "/home/Project/Metabolome/HFQ/HFQ_Final-1/P_Dis/ID_Metabolites/FDN.ID.txt", "/home/Project/Metabolome/HFQ/HFQ_Final-1/P_Dis/bin/FDN-grouping.info" ,"/home/Project/Metabolome/HFQ/HFQ_Final-1/P_Dis/4.10_ID_KEGG" )
setwd(arg[3])


library(MetaboAnalystR)
###
id <- read.table(arg[1],header=T,sep='\t')
n.name <- strsplit(basename(arg[1]),".",fixed=T)[[1]][1]

mSet<-InitDataObjects("conc", "pathora", FALSE)
cmpd.vec <- as.character(id[,2])

mSet <- Setup.MapData(mSet, cmpd.vec)
mSet <- CrossReferencing(mSet, "hmdb")
mSet <- SetKEGG.PathLib(mSet, "hsa")    #human
#

nm.map <- GetFinalNameMap(mSet)
write.csv(nm.map,file = paste0(n.name,".kegg.csv"),row.names=F,col.names=T,quote=F)
valid.inx <- !(is.na(nm.map$kegg) | duplicated(nm.map$kegg))
ora.vec <- nm.map$kegg[valid.inx]
q.size <- length(ora.vec)
#if (is.na(ora.vec) || q.size == 0) {
#    AddErrMsg(mSet, "No valid KEGG compounds found!")
#    return(0)
#}
library(KEGGgraph)
current.mset <- metpa$mset.list
uniq.count <- metpa$uniq.count
hits <- lapply(current.mset, function(x) {
    x[x %in% ora.vec]
})
hit.num <- unlist(lapply(hits, function(x) {
    length(x)
}), use.names = FALSE)
set.size <- length(current.mset)
set.num <- unlist(lapply(current.mset, length), use.names = FALSE)
res.mat <- matrix(0, nrow = set.size, ncol = 8)
rownames(res.mat) <- names(current.mset)
colnames(res.mat) <- c("Total", "Expected", "Hits", "Raw p",
"-log(p)", "Holm adjust", "FDR", "Impact")
imp.list <- metpa$rbc      #relative betweenness centrality
#imp.list <- metpa$dgr      #out degree centrality
res.mat[, 1] <- set.num
res.mat[, 2] <- q.size * (set.num/uniq.count)
res.mat[, 3] <- hit.num
res.mat[, 4] <- GetFisherPvalue(hit.num, q.size, set.num,uniq.count)  #Fishers' exact test
#res.mat[, 4] <- phyper(hit.num - 1, set.num, uniq.count-set.num, q.size, lower.tail = F)   #Hypergeometric test
res.mat[, 5] <- -log(res.mat[, 4])
res.mat[, 6] <- p.adjust(res.mat[, 4], "holm")
res.mat[, 7] <- p.adjust(res.mat[, 4], "fdr")
res.mat[, 8] <- mapply(function(x, y) {
    sum(x[y])
}, imp.list, hits)
res.mat <- res.mat[hit.num > 0, , drop = FALSE]
if (nrow(res.mat) > 1) {
    ord.inx <- order(res.mat[, 4], res.mat[, 8])
    res.mat <- res.mat[ord.inx, ]
}
mSet$analSet$ora.mat <- signif(res.mat, 5)


mSet$analSet$ora.hits <- hits
save.mat <- mSet$analSet$ora.mat
rownames(save.mat) <- GetORA.pathNames(mSet)
write.csv(save.mat, file = paste0(n.name,".MetaboAnalystR_pathway_results.csv"))

Result <- cbind(rownames(mSet$analSet$ora.mat),GetORA.pathNames(mSet))
colnames(Result) <- c("ID1","Pathway")
rr <- matrix(NA,ncol=2,nrow=0)
for(i in 1:nrow(Result)){
	r.dd <- hits[[as.character(Result[i,1])]]
	r1 <- paste0(names(r.dd),collapse="/")
	r2 <- paste0(r.dd,collapse="/")
	r3 <- cbind(r1,r2)
	rr <- rbind(rr,r3)
}
colnames(rr) <- c("Metabolitics","ID2")
Result <- cbind(Result,rr)
write.table(Result,file = paste0(n.name,".Result.Kegg.txt"),quote=F,sep="\t",row.names=F,col.names=T)


#######
pdf(paste0(n.name,".MetaboAnalystR_pathway.pdf"),width=10,height=8);
layout(mat=matrix(c(1,2,1,3),nrow=2,byrow=T),w=c(8,1),h=c(1,1))
logp <- as.numeric(save.mat[,'-log(p)'])
path.impact <- as.numeric(save.mat[,'Impact'])
cPal <- colorRampPalette(c("LemonChiffon","Firebrick"))
mycolors <- paste0(cPal(length(logp))[as.numeric(cut(logp,breaks = length(logp)))],'')

gg <- read.table(arg[2],header=T,row.names=1)
m.lab <- unique(as.character(gg[,1]))

par(mgp=c(3,1,0),mar=c(6,6,3.5,6),cex.axis=2,lwd=1,font=1,adj=0.5,cex.lab=2,cex.main=2)

plot(0,0,xlab="Pathway Impact",ylab="-log(p)",axes=F,main=paste0(m.lab[1]," vs ",m.lab[2]),xlim=c(0,max(path.impact)*1.1),ylim=c(min(logp)*0.9,round(max(logp))*1.1),col="white",xpd=T)
axis(side=1,at=signif(seq(0,max(path.impact),max(path.impact)/3),2))
axis(side=2,at=round(seq(round(min(logp)),round(max(logp)),length.out=4)),las=2)

p.cex <- path.impact
un.p <- unique(p.cex)
p.pro <- round(quantile(un.p,c(0.25,0.5,0.75)),3)
p.t <- p.cex
p.cex[which(p.t <= p.pro[1])] <-1
p.cex[which(p.t <= p.pro[2] & p.t > p.pro[1])] <-2
p.cex[which(p.t <= p.pro[3] & p.t > p.pro[2])] <-3
p.cex[which(p.t > p.pro[3])] <-4


points(path.impact,logp,col=mycolors,pch=16,xpd=T,cex=1.5*p.cex,xpd=T)
points(path.impact,logp,col='black',pch=1,lwd=0.00001,xpd=T,cex=1.5*p.cex,xpd=T)
#text(path.impact[path.impact >0],logp[path.impact >0],labels=rownames(save.mat)[path.impact >0],xpd=T,cex=0.8,font=2,xpd=T)
ord <- order(logp,decreasing=T)
#ord <- ord[1:5]
xx <- path.impact[ord]
yy <- logp[ord]
ll <- rownames(save.mat)[ord]
#nn <- floor(length(xx)/2)]
nn <-5
text(xx[1:nn],yy[1:nn],labels=ll[1:nn],xpd=T,cex=0.8,font=1,xpd=T)

ll1 <- "Glycerophospholipid metabolism"
if(!(ll1 %in% ll)){
	ord <- save.mat[ll1,]
	xx <- ord["Impact"]
	yy <- ord["-log(p)"]
	text(xx,yy,labels=ll,xpd=T,cex=0.8,font=1,xpd=T)
}


#abline(v=signif(seq(0,max(path.impact),length.out=5),2),lty=2,col='blue')
#abline(h=round(seq(round(min(logp)),round(max(logp)),length.out=5)),lty=2,col='blue')


##	image
par(mar=c(2,0,4,4))
dd <- res.mat[,"Raw p"]
d.or <- order(dd,decreasing=T)
dd <-dd[d.or]
image(t(as.matrix(seq(max(dd),min(dd),length.out=length(mycolors)))),axes=F,col=mycolors,xlim=c(0,1))
axis(4,at=seq(0,1,length.out=3),labels=round(c(seq(max(dd),min(dd),length.out=3)),1),tick=F,line=0,cex.axis=1.2,las=1)
#mtext(expression('log'[10]*'(Abundance)'),side=3,line=0,cex=1.2)
mtext('P.Value',side=3,line=0,cex=1.2)	

##	legend
par(mar=c(6,0,9,4),cex.main=2,cex.main=1.2)
p.p <- unique(p.cex)
ord <- order(p.p,decreasing=F)
p.p <- p.p[ord]
plot(rep(0.5,length(p.p)),p.p,col="grey",xlim=c(0,1),xlab="",ylab="",axes=F,cex=1.5*p.p,pch=16,xpd=T)
p.ll <- c(paste0("< ",p.pro[1]),paste0("< ",p.pro[2]),paste0("< ",p.pro[3]),paste0("> ",p.pro[3]))
#text(rep(1,length(p.p)),p.p,labels=p.ll,xpd=T,cex=0.8,font=1,xpd=T)
#legend("left",fill=unique(PopCol),legend=unique(g.c[,1]),horiz=T,xpd=T,bty="n")
mtext('Pathway Impact',side=3,line=1,cex=1.2)	
mtext(text=p.ll,side=4,cex=1.2,at=p.p,las=1)
dev.off()


