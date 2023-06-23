#Analysis StepI: Classify the abundance data 
#Current working taxanomic level is Genus

Data=$(dirname $PWD)

[ ! -d $Data/Analysis ] && mkdir $Data/Analysis;
[ ! -d $Data/Analysis/3_PCA_cluster ] && mkdir $Data/Analysis/3_PCA_cluster;

cp $Data/Analysis/1_Abundance/*.Genus.abundance.txt $Data/Analysis/3_PCA_cluster/;

for i in $(ls $Data/Analysis/3_PCA_cluster/*.Genus.abundance.txt)
do
	wd=$(dirname $i)
	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	echo $c1
	Rscript $Data/bin/A3.1_Plot_PCA_cluster.r $wd $i 0.001 0 0 $Data/bin/$c1-color.txt;
	#参数：工作空间；数据；
	## 部分数据画不出
done

#bash A3.0_Extract_Group_abundance.sh
#bash A3.2_Group_figures.sh
