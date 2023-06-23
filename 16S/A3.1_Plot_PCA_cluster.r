# Reference:http://enterotype.embl.de/enterotypes.html
# Please make sure packages, including ade4, cluster, clusterSim, gplots, and MASS are already installed.
# Every column in source data stands for an individual.
#
# Suggested noise.removal percent is 0.01.
# There will be 6 output pdf files for abc.txt, including abc_cluster.pdf, abc_cluster_cluster.txt, abc_cluster_eva.pdf ,abc_top_three_box.pdfand abc_top_three_bar.pdf. 
print("This script should be used as cluster_plot_v3.0.R working_dir(complete path) file_name noise.removal_percent removal_number k-number")
print("Set k as 0 for automatically choosing best k number.")
print("If this doesn't work, check comments.")
args <- commandArgs(TRUE)
if(length(args)<5 | length(args)>7){
    cat("Argument: working_dir file_name noise.removal_percent removal_number k-number [ Distmethod Validmethod ]\n")
    cat("Distmethod: rJSD (as default) JSD bray\nValidmethod: CH (as default) PS (recommend) SI\n")
    quit('no')
}

print(args)
#args <-c("/home/Project/16S/JZM/Analysis/3_PCA_cluster"   , "/home/Project/16S/JZM/Analysis/3_PCA_cluster/Total.Genus.abundance.txt", "0.001" , "0"   ,"0"  )



file_dir <- as.character(args[1])   # read working directory 
file_name <- as.character(args[2])   # read file name
percent <- as.numeric(args[3])   # read noise.reoval percent 
rm_n <- as.numeric(args[4])    # set remove number of top classes
k <- as.numeric(args[5])  # set k number
grp=read.table(args[6],sep='\t',header=T,row.names=1,comment.char = "")

Krange = 10
Distmethod = 'rJSD'
Validmethod = 'CH'
if(length(args)==7){Distmethod = args[7]}
if(length(args)==8){Validmethod = args[8]}

setwd(file_dir)   #set working directory

# load required packages
suppressMessages(require(MASS))
suppressMessages(require(cluster))
suppressMessages(require(clusterSim))
suppressMessages(require(ade4))
suppressMessages(require(fpc))
suppressMessages(require(vegan))

KLD <- function(x,y) {sum(x *log(x/y))}
JSD<- function(x,y) {(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))}
rJSD<- function(x,y) {sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))}
# generate noise.removal function
noise.removal <- function(dataframe, percent=0.01, top=NULL){
        dataframe->Matrix
        bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
        Matrix_1 <- Matrix[bigones,]
        print(paste("noise.removal percent is", percent, sep = " "))
        return(Matrix_1)
}
# genertae JSD function
dist <- function(inMatrix, method = 'rJSD', pseudocount=0.000001, ...){
    matrixColSize <- length(colnames(inMatrix))
    matrixRowSize <- length(rownames(inMatrix))
    colnames <- colnames(inMatrix)
    inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

    if(method == 'rJSD'){
        resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
        for(i in 1:matrixColSize) {
            for(j in 1:matrixColSize) {
                resultsMatrix[i,j]=rJSD(as.vector(inMatrix[,i]),as.vector(inMatrix[,j]))
            }
        }
        colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
        as.dist(resultsMatrix)->resultsMatrix
        attr(resultsMatrix, "method") <- "rJSD"
    }else if(method == 'JSD'){
        resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
        for(i in 1:matrixColSize) {
            for(j in 1:matrixColSize) {
                resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),as.vector(inMatrix[,j]))
            }
        }
        colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
        as.dist(resultsMatrix)->resultsMatrix
        attr(resultsMatrix, "method") <- "JSD"
    }else if(method == 'bray'){
        require(vegan)
        resultsMatrix <- vegdist(t(inMatrix))   # column indicates taxa
    }else{
        cat("Calculate distance with 'rJSD', 'JSD',or 'bray'.\n")
        quit('no')
    }
    return(resultsMatrix)
}
pam.clustering <- function(x,k) { # x is a distance matrix and k the number of clusters
        cluster <- as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
        return(cluster)
}
best_cluster <- function(nclu) {
        best_clu <- 2
        for (i in 2: 20)
                if (nclu[i] > nclu[best_clu]) {best_clu <- i}
        return(best_clu)
}
max_name <- function(clu) {
        max_col <- 1
        for (i in 1:(ncol(clu)-1)){
                if (clu[nrow(clu), i] > clu[nrow(clu), max_col]) max_col <- i 
        }
        return(colnames(clu)[max_col])
}
clu_bar_three <- function(clu){
        clu_bar <- clu[nrow(clu),][colnames(clu) == rank_name(clu)[1]]
        for (m in 2:3) {
                clu_bar <- cbind(clu_bar, clu[nrow(clu),][colnames(clu) == rank_name(clu)[m]])        
        }
        return(clu_bar)
}
rank_name <- function(clu) {
        rank1 <- max_name(clu)
        clu_tmp <- clu[colnames(clu) != rank1]
        rank2 <- max_name(clu_tmp)
        rank3 <- max_name(clu_tmp[colnames(clu_tmp) != rank2])
        return(c(rank1, rank2, rank3))
}
clu_box_three <- function(clu){
        clu_box <- clu[colnames(clu) == rank_name(clu)[1]]
        for (m in 2:3) {
                clu_box <- cbind(clu_box, clu[colnames(clu) == rank_name(clu)[m]])        
        }
        return(clu_box)
}
clu_label <- function(clu) {
#        clu_label <- paste(rank_name(clu)[1],"\n", rank_name(clu)[2], sep="")  ## first two 
        clu_label <- rank_name(clu)[1]
        return(clu_label)
}
rm_top_class <- function(data) {
        data <- t(data)
        data_tmp <- rbind(data, apply(data, 2, mean))
        max_class <- 1
        for (i in 1:(ncol(data_tmp))){
                if (data_tmp[nrow(data_tmp), i] > data_tmp[nrow(data_tmp), max_class]) max_class <- i 
        }
        return(t(data[, -max_class]))
}

OptimizeK <- function(data,method='CH',Kmin=1,Kmax=Krange, file_name){
    data.dist <- dist(data,method=Distmethod)
    obs.silhouette <- c()
    nclusters <- c()
    si.optimalk <- 1
    CH.optimalk <- 1
    for(k in Kmin:Kmax){
        data.cluster_temp = pam.clustering(data.dist,k)
        if(k > 1){
            nclusters[k] = index.G1(t(data), data.cluster_temp, d = data.dist, centrotypes = "medoids")
            if(!is.na(nclusters[k])){
                if(is.na(nclusters[CH.optimalk])){
                    CH.optimalk = k
                }else{
                    CH.optimalk = ifelse(nclusters[k]>=nclusters[CH.optimalk],k,CH.optimalk)
                }
            }
            obs.silhouette[k] <- mean(silhouette(data.cluster_temp, data.dist)[,3])
            if(!is.na(obs.silhouette[k])){
                if(is.na(obs.silhouette[si.optimalk])){
                    si.optimalk = k
                }else{
                    si.optimalk = ifelse(obs.silhouette[k]>=obs.silhouette[si.optimalk],k,si.optimalk)
                }
            }
        }else{
            nclusters[k] = NA
            obs.silhouette[k] = NA
        }
    }
    ps=prediction.strength(data.dist,Gmin=ifelse(2>Kmin, 2, Kmin),Gmax=Krange,M=50,clustermethod=pamkCBI)

    cat("Distance method:\t",Distmethod,"\n")
    cat("PS validation:\nmean.PS:\t",ps$mean.pred,"\ncutoff:\t",ps$cutoff,"\noptimalK:\t",ps$optimalk,"\n")
    cat("SI validation:\nSI values:\t",obs.silhouette,"\noptimalK:\t",si.optimalk,"\n")
    cat("CH validation:\nCH values:\t",nclusters,"\noptimalK:\t",CH.optimalk,"\n")

    pdf(paste(file_name, "_cluster_eva.pdf", sep=""))
    plot(nclusters, type="h", xlab="j clusters", ylab="CH index")
    dev.off()

    if(Validmethod == 'PS'){
        optimalK = ps$optimalk
    }else if(Validmethod == 'SI'){
        optimalK = si.optimalk
    }else{
        optimalK = CH.optimalk
    }
    return(optimalK)
}

cluster_plot <- function(data, file_name_rt, k , percent,grp) {

        data.dist <- dist(data,method=Distmethod)   # generate distance matrix

        if(k == 0){
            k <- OptimizeK(data,method = Validmethod, Kmin=1,Kmax=Krange, file_name= file_name_rt)
        }else{
            OptimizeK(data,method = Validmethod, 1, Krange, file_name= file_name_rt)
        }
        
        # use pam to obtain cluster result
        data.cluster <- pam.clustering(data.dist, k)
        obs.silhouette <- mean(silhouette(data.cluster, data.dist)[,3])
        
        write.table(data.cluster, paste(file_name_rt, "_cluster.txt", sep= ""), row.names =T)
        
        
        # generate coordinates for each individual
        #obs.pca <- dudi.pca(data.frame(t(data)), scannf=F, nf=10)
        #obs.bet <- bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1)    
        
        obs.pcoa=dudi.pco(data.dist, scannf=F, nf=k)
        
        
        # assign pch parameter to different races
        data.sample <- obs.pcoa$li
        n <- ncol(data.sample)
		gg0=unique(grp)
		gg=as.character(gg0[,1])
		grp_col=as.character(gg0[,2])
		names(grp_col)=gg
		pch=rep(1:25)
		names(pch)=gg
		dd.s=rownames(data.sample)
		dd.s0=as.character(grp[dd.s,1])
		dd.s=as.numeric(pch[dd.s0])
		data.sample[, n+1]=dd.s
		dd.c=grp_col[dd.s0]
        colnames(data.sample)[ncol(data.sample)] <- "race"
        
        
        data.clu <- data.frame(t(data), data.cluster)
        colnames(data.clu)[ncol(data.clu)] <- "clu_type"
        
        write.table(data.frame(rownames(data.clu), data.cluster), paste(file_name_rt, "_cluster.txt", sep= ""))
                
        clu_names <- NULL
        clu_labels <- NULL
        top_three <- NULL
        clu_bat_names <- NULL
        clu_bot_names <- NULL
        
        for (j in 1:k) {  
                clu_names <- rbind(clu_names, paste("cluster", j, sep = ""))
                assign(clu_names[j], NULL)
                clu_bat_names <- rbind(clu_bat_names, paste("clu", j, "_bat", sep=""))
                assign(clu_bat_names[j], NULL)
                clu_bot_names <- rbind(clu_bot_names, paste("clu", j, "_bot", sep=""))
                assign(clu_bot_names[j], NULL)
                for (i in 1:nrow(data.clu)) {
                        if (data.clu[i, ncol(data.clu)] == j) assign(clu_names[j], rbind(get(clu_names[j]), data.clu[i, ]))         
                }
                assign(clu_names[j], rbind(get(clu_names[j]), apply(get(clu_names[j]), 2, mean)))
                clu_labels <- rbind(clu_labels, clu_label(get(clu_names[j])))
                top_three <- rbind(top_three, rank_name(get(clu_names[j])))
                assign(clu_bat_names[j], clu_bar_three(get(clu_names[j])))
                assign(clu_bot_names[j], clu_box_three(get(clu_names[j])))
                
        }
        
        pdf(paste(file_name_rt, "_top_three_bar.pdf", sep = ""), width = 18, height = 5)
        layout(matrix(1:k, 1, k, byrow = T))
        for (j in 1:k) {
                barplot(as.matrix(get(clu_bat_names[j])), main = paste("Top 3 classes in ", clu_names[j], sep=""), beside = T, cex.names = 1.3, xlab = "Classes", ylab = "Relative abundance", col = c("darkslategray1", "gold", "dodgerblue"))  
        }
        dev.off()
        
        pdf(paste(file_name_rt, "_top_three_box.pdf", sep = ""), width = 15, height = 5)
        layout(matrix(1:k, 1, k, byrow = T))
        for (j in 1:k) {
                boxplot(as.matrix(get(clu_bot_names[j])), main = paste("Top 3 classes in ", clu_names[j], sep=""), beside = T, cex.names = 1.3, xlab = "Classes", ylab = "Relative abundance", col = c("darkslategray1", "gold", "dodgerblue"))  
        }
        
        dev.off()
        
        
        pdf(paste(file_name_rt, "_cluster.pdf", sep = ""))
        eig_pc1_contr = round(obs.pcoa$eig[1]/sum(obs.pcoa$eig)*100,2)
        eig_pc2_contr = round(obs.pcoa$eig[2]/sum(obs.pcoa$eig)*100,2)
 #       plot(obs.pcoa$li[,1:2], type = "n", frame.plot = T, xlab = paste("PC1(",eig_pc1_contr,"%)",sep=''), ylab = paste("PC2(",eig_pc2_contr,"%)",sep=''))
		  plot(obs.pcoa$li[,1:2], frame.plot = T, xlab = paste("PC1(",eig_pc1_contr,"%)",sep=''), ylab = paste("PC2(",eig_pc2_contr,"%)",sep=''),pch=dd.s,col=dd.c)
		

        abline(h = 0, v = 0)
#        s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, addaxes = T, pch = data.sample$race, label = clu_labels, clabel = 1, axesell =F, col = clu_colors, cpoint = 1.2, add.plot = T)
        title(paste(file_name_rt, "_cluster", sep = ""))
        legend("topleft", title = "Races Legends", pch = pch[gg], col=grp_col[gg],legend = gg,, cex = 1)
#        legend("bottomright", title = "Cluster Legend", legend = clu_names, cex = 1.3, col= clu_colors, pch = 18)
        dev.off()
        
        top_three <- data.frame(top_three)
        rownames(top_three) <- clu_names
        colnames(top_three) <- c("Rank1", "Rank2", "Rank3")
        write.table(top_three, paste(file_name_rt, "_top_three.txt", sep= ""))
        print(paste("Script runs successfully", " for ", file_name_rt, ".", sep = ""))      
}
rm_cluster_plot <- function(file_name, k , percent, rm_n,grp) {
        data <- read.table(file_name, header=T, row.names=1, dec=".",check.names = FALSE)   # read abundance matrix
        file_name_rt <- strsplit(file_name, ".txt")
        data <- noise.removal(data, percent)  # genertae denoized data
        for (i in 1:(rm_n+1)) {
                cluster_plot(data, file_name_rt, k , percent,grp)
                data <- rm_top_class(data)
                file_name_rt <- paste(file_name_rt, "_rm", i, sep = "")       
        }     
}

rm_cluster_plot(file_name, k , percent, rm_n,grp)
