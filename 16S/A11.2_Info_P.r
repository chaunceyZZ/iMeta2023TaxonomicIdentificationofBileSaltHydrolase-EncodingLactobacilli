arg1 <- commandArgs(T)
if(length(arg1) < 1){
	cat("Argument: Out_Dir \n")
	quit('no')
}

print(arg1)
###########	测试
#arg1 <- c( "/home/wmj/视频/16s/qiime1/TT/Analysis/15_RDA_and_Heatmap/Total", "T1Control-T1Classical"     ,"T1Control-T1T1H" , "T1Control-T1T1L","0.05"  , "3", "/home/wmj/视频/16s/qiime1/TT/bin/Total-color.txt"  )

setwd(arg1[1])
len.arg1 <- length(arg1)-1
arg.tt <- arg1[2:(length(arg1)-3)]
arg <- arg.tt
###########	自定义函数区
### 去除NaN
F_del_na <- function(data) {
	t1 <- sum(is.na(data))
	return(t1)
}

### P值符合要求
P <- as.numeric(arg1[len.arg1-1])
len.fp <- as.numeric(arg1[len.arg1])
F_p_save <- function(data) {
	nn <- which(data <= P)
	if(length(nn) == 0) return(0)    #### 所有P >0.05
	if(length(nn) >= len.fp) {       #### 所有P <=0.05
		return(len.fp)
	} else {	 #### 介于两者直接
		return(0.5)
	}
}


#####################################################
grpInfo = read.table(file=arg1[length(arg1)],header=T,sep="\t",check.names = FALSE, comment.char = "")
rownames(grpInfo)=as.character(grpInfo[,1])
pvalue <- read.csv(file="02_p_value.csv",header=T,row.names=1)
gut <- read.csv(file="02_all_gut.csv",header=T,row.names=1)
ii=intersect(rownames(grpInfo),colnames(gut))
if(length(ii)!=nrow(grpInfo)){
	print('Data Lost')
	ii1=setdiff(rownames(grpInfo),ii)
	print(ii1)
}
grpInfo=grpInfo[ii,]
gut=gut[,ii]

gut_col <- colnames(gut)

grp <- length(arg)

########## 比较的P值
p_col <- colnames(pvalue)
cc <- c()
for (i in 1:grp) {
	group <- unlist(strsplit(arg[i],split="-"))
	m1 <- paste(group[1],".*",group[2],sep="")
	c1 <- grep(m1,p_col,value=F,fixed = F,ignore.case = F)
	cc <- c(cc,c1)
}

#######	挑出data	
## del na
if (grp != 1){
	len_na <- apply(pvalue[,cc],1,F_del_na)

} else {
	len_na <- which(is.na(pvalue[,cc]))
}
if(length(len_na)>0){
	p1 <- as.matrix(pvalue[-which(len_na>0),cc])
	all_p1 <- as.matrix(pvalue[-which(len_na>0),])
	g1 <- gut[-which(len_na>0),]
} else {
	p1 <- as.matrix(pvalue[,cc])
	all_p1 <- as.matrix(pvalue[,])
	g1 <- gut
}

######## 挑选p
len_p <- apply(p1,1,F_p_save)
p2 <- as.matrix(p1[which(len_p > 0),])
all_p2 <- as.matrix(all_p1[which(len_p > 0),])
g2 <- g1[which(len_p >0),]

rownames(p2) <- rownames(g2)
rownames(all_p2) <- rownames(g2)
colnames(p2) <- colnames(p1)
colnames(all_p2) <- colnames(pvalue)

write.csv(p2,file="03_p_value.csv",row.names=T)
write.csv(all_p2,file="03_all_p_value.csv",row.names=T)
write.csv(g2,file="03_all_gut.csv",row.names=T)






