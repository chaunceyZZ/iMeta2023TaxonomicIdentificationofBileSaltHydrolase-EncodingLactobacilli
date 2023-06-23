arg <- commandArgs(T)
if(length(arg) != 1){
	cat("Argument: related.abundance.txt\n")
	quit('no')
}
filepath <- strsplit(arg[1],'/',fixed=T)[[1]]
filename <- filepath[length(filepath)]
level <- strsplit(filename,'.',fixed=T)[[1]][2]

Abundance <- read.table(arg[1],header=T)
HAN.abundance <- as.matrix(data.frame(Abundance[,grep("^H\\.",names(Abundance))]))
KZK.abundance <- as.matrix(data.frame(Abundance[,grep("^K\\.",names(Abundance))]))
UIG.abundance <- as.matrix(data.frame(Abundance[,grep("^U\\.",names(Abundance))]))

rownames(HAN.abundance)=Abundance[,1]
rownames(KZK.abundance)=Abundance[,1]
rownames(UIG.abundance)=Abundance[,1]

#####Calculate the total abundance at each individual and boxplot the distribution
HAN.indi.sum.abundance=colSums(HAN.abundance);
KZK.indi.sum.abundance=colSums(KZK.abundance);
UIG.indi.sum.abundance=colSums(UIG.abundance);
pdf(paste("Own.",as.character(level),".indi.sum.abundance.pdf",sep=''),width=5,height=5);
boxplot(HAN.indi.sum.abundance,KZK.indi.sum.abundance,UIG.indi.sum.abundance,border=c("red","blue","green"),names=c("HAN","KZK","UIG"));
dev.off()
Indi.sum.abundance=c(HAN.indi.sum.abundance,KZK.indi.sum.abundance,UIG.indi.sum.abundance);
write.csv(Indi.sum.abundance,file=paste("Own.",as.character(level),".indi.sum.abundance.txt",sep=''));

#######Calculate the abundance of each type and boxplot the distribution
pdf(paste("Own.",as.character(level),".type.pop.abundance.pdf",sep=''),width=10,height=7);
#par(mfrow=c(1,1),mar=c(30,3,3,3));
Type.abundance=as.matrix(data.frame(Abundance[,2:ncol(Abundance)],sum=rowSums(Abundance[,2:ncol(Abundance)])));
rownames(Type.abundance)=Abundance[,1];
Type.abundance=Type.abundance[order(Type.abundance[,ncol(Type.abundance)],decreasing=T),];	# sorted by the total abundance
Type.abundance=Type.abundance[,-ncol(Type.abundance)];
HAN.abundance=Type.abundance[,grep("^H\\.",names(Abundance))-1];
KZK.abundance=Type.abundance[,grep("^K\\.",names(Abundance))-1];
UIG.abundance=Type.abundance[,grep("^U\\.",names(Abundance))-1];

Type.abundance.top = Type.abundance;
if(nrow(Type.abundance.top)>15){															# get top 15 types
	Type.abundance.top=Type.abundance.top[-c(16:nrow(Type.abundance.top)),];
}

HAN.abundance.top=Type.abundance.top[,grep("^H\\.",names(Abundance))-1];
KZK.abundance.top=Type.abundance.top[,grep("^K\\.",names(Abundance))-1];
UIG.abundance.top=Type.abundance.top[,grep("^U\\.",names(Abundance))-1];
Type.group=matrix(,nrow=nrow(Type.abundance.top),ncol=ncol(Type.abundance.top));
Type.color=matrix(,nrow=nrow(Type.abundance.top),ncol=ncol(Type.abundance.top));
Type.name=rownames(Type.abundance.top);
for(m in 1:nrow(Type.abundance.top)){
	Type.group[m,]=c(rep(3*(m-1)+1,ncol(HAN.abundance.top)),rep(3*(m-1)+2,ncol(KZK.abundance.top)),rep(3*(m-1)+3,ncol(UIG.abundance.top)));
	Type.color[m,]=c(rep("red",ncol(HAN.abundance.top)),rep("blue",ncol(KZK.abundance.top)),rep("green",ncol(UIG.abundance.top)));
}

#boxplot(Type.abundance.top~Type.group,border=c("red","blue","green"),xaxt="n");
#rect(seq(0.5,max(Type.group),6),-0.1,seq(3.5,max(Type.group)+0.5,6),1,col='grey90',lty=0)
#boxplot(Type.abundance.top~Type.group,border=c("red","blue","green"),xaxt="n",add=T);
#axis(1,labels=Type.name,at=seq(2,44,3),las=2); 

####perform wilcox U test to compare the abundance of each type among three populations####
Results=matrix(,nrow=nrow(Type.abundance),ncol=6);
rownames(Results)=rownames(Type.abundance);
colnames(Results)=c("HAN.AVG","KZK.AVG","UIG.AVG","HAN.KZK.U.pBH","HAN.UIG.U.pBH","KZK.UIG.U.pBH");
Results[,1]=apply(HAN.abundance,1,mean);
Results[,2]=apply(KZK.abundance,1,mean);
Results[,3]=apply(UIG.abundance,1,mean);
Type.significant=rownames(Type.abundance.top);
for(i in 1:nrow(Type.abundance)){
	Results[i,4]=wilcox.test(as.vector(as.matrix(HAN.abundance[i,])),as.vector(as.matrix(KZK.abundance[i,])))$p.value;
	Results[i,5]=wilcox.test(as.vector(as.matrix(HAN.abundance[i,])),as.vector(as.matrix(UIG.abundance[i,])))$p.value;
	Results[i,6]=wilcox.test(as.vector(as.matrix(UIG.abundance[i,])),as.vector(as.matrix(KZK.abundance[i,])))$p.value;
}
Results[,4] = p.adjust(Results[,4],'BH');
Results[,5] = p.adjust(Results[,5],'BH');
Results[,6] = p.adjust(Results[,6],'BH');

write.csv(Results,file=paste("Own.",as.character(level),".pop.abundance.utest.BH.txt",sep=''));

top_top_n = 0;
for(r in 1:nrow(Type.abundance.top)){
if(sum(Results[r,1:3])>=0.05){top_top_n = r}
else{break}
}
rect_bot = seq(0.5,max(Type.group),3*2)
rect_top = seq(3+0.5,max(Type.group)+0.5,3*2)
par(mfrow=c(1,1),mar=c(20,3,4,3));
boxplot(Type.abundance.top[1:top_top_n,]~Type.group[1:top_top_n,],col=c("red","blue","green"),xaxt="n",xlim=c(1,max(Type.group)),outline=F);
rect(rect_bot[1:ceiling(top_top_n/2)],-0.1,rect_top[1:ceiling(top_top_n/2)],1,col='grey90',lty=0)
boxplot(Type.abundance.top[1:top_top_n,]~Type.group[1:top_top_n,],col=c("red","blue","green"),xaxt="n",add=T,outline=F);
par(new=T)
boxplot(Type.abundance.top[(top_top_n+1):nrow(Type.abundance.top),]~Type.group[(top_top_n+1):nrow(Type.abundance.top),],col=c("red","blue","green"),xaxt="n",xlim=c(1,max(Type.group)),at=min(Type.group[(top_top_n+1):nrow(Type.abundance.top),]):max(Type.group[(top_top_n+1):nrow(Type.abundance.top),]),yaxt='n',outline=F);
rect(rect_bot[(ceiling(top_top_n/2)+1):length(rect_bot)],-0.1,rect_top[(ceiling(top_top_n/2)+1):length(rect_top )],1,col='grey90',lty=0)
boxplot(Type.abundance.top[(top_top_n+1):nrow(Type.abundance.top),]~Type.group[(top_top_n+1):nrow(Type.abundance.top),],col=c("red","blue","green"),add=T,at=min(Type.group[(top_top_n+1):nrow(Type.abundance.top),]):max(Type.group[(top_top_n+1):nrow(Type.abundance.top),]),yaxt='n',xaxt='n',outline=F);
abline(v=max(Type.group[1:top_top_n,])+0.5,lty=2,lwd=2)
axis(4)
axis(1,labels=Type.name,at=seq(2,44,3),las=2)
#mtext(side=1,at=seq(2,44,3),text=Type.name,srt=270)
mtext(side=3,line=1,at=c(mean(rect_bot[1:ceiling(top_top_n/2)]),mean(rect_bot[(ceiling(top_top_n/2)+1):length(rect_bot)])),text=c(expression(sum(AGV_RA,pops)>='5%'),expression(sum(AGV_RA,pops)<'5%')),cex=0.8)

for(i in 1:nrow(Type.abundance.top)){
	Type.significant[i]="   ";
	if(min(Results[i,4],Results[i,5],Results[i,6])<0.05){
		Type.significant[i]=" * ";
	}
	if(min(Results[i,4],Results[i,5],Results[i,6])<0.01){
		Type.significant[i]=" ** ";
	}
	if(min(Results[i,4],Results[i,5],Results[i,6])<0.001){
		Type.significant[i]="***";
	}
}
axis(3,labels=Type.significant,at=seq(2,44,3),tick = F,line = -1.2);
dev.off()

