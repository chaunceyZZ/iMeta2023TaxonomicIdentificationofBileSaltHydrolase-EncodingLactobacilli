#Analysis StepI: Classify the abundance data 
#Current working taxanomic level is Genus

Data=$(dirname $PWD)

[ ! -d $Data/Analysis/5_Sharing_OTUs ] && mkdir $Data/Analysis/5_Sharing_OTUs;

cp $Data/Analysis/4_common_OTUs/*/*.otu.counts.txt $Data/Analysis/5_Sharing_OTUs/;

for i in $(ls $Data/Analysis/5_Sharing_OTUs/*.otu.counts.txt)
do
	Out_dir=$(dirname $i)
	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	Rscript $Data/bin/A5.1_Plot_OTUs_Venn.r $i $Data/bin/$c1-color.txt $Out_dir 2;

	#参数：Data_File Group_File Out_Dir Group_column
done


