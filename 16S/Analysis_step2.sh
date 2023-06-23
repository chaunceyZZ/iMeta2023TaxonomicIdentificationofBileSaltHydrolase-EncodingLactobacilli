#Analysis StepI:Convert the raw Qiime abundance file to the tab-format files!
#Current working taxanomic level is Genus

Data=$(dirname $PWD)

[ ! -d $Data/Analysis ] && mkdir $Data/Analysis;
[ ! -d $Data/Analysis/2_Diversity ] && mkdir $Data/Analysis/2_Diversity;

cp $Data/Analysis/1_Abundance/*.Genus.abundance.txt $Data/Analysis/2_Diversity/;
for i in $(ls $Data/Analysis/2_Diversity/*.Genus.abundance.txt)
do
	outdir=$(dirname $i)
	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	Rscript $Data/bin/A2.1_Plot_alpha_bar.r $i $Data/bin/$c1-color.txt $outdir ;
	Rscript $Data/bin/A2.2_Plot_beta_bar.r $i $Data/bin/$c1-color.txt $outdir ;
	#参数：数据文件；组别信息；工作目录
done

cd $Data/bin
#bash A2.0_Extract_Group_abundance.sh
#bash A2.3_Group_figures.sh
