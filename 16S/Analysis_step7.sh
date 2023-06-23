#Analysis StepI: Classify the abundance data 
#Current working taxanomic level is Genus

Data=$(dirname $PWD)


###	7_OTUs_alpha
for i in $(ls $Data/14_Statistics/alpha/*/*final.otu.*.txt)
do
	dir1=${i#*alpha}
	dir2=`echo "$dir1" | cut -d / -f 2`
	[ ! -d $Data/Analysis/7_OTUs_alpha/$dir2 ] && mkdir -p $Data/Analysis/7_OTUs_alpha/$dir2;
	cp $i $Data/Analysis/7_OTUs_alpha/$dir2
	Rscript $Data/bin/A7.1_Plot_OTUs_alpha_bar_modify.r $i $Data/bin/$dir2-color.txt $Data/Analysis/7_OTUs_alpha/$dir2 2;	
done


###	8_OTUs_beta 


for i in $(ls $Data/14_Statistics/beta/*/*.otu_table*.txt);
do
	dir1=${i#*beta}
	dir2=`echo "$dir1" | cut -d / -f 2`
	[ ! -d $Data/Analysis/8_OTUs_beta/$dir2 ] && mkdir -p $Data/Analysis/8_OTUs_beta/$dir2;
	cp $i $Data/Analysis/8_OTUs_beta/$dir2;
	Rscript $Data/bin/A7.3_Plot_OTUs_beta_mds.r  $i $Data/bin/$dir2-color.txt $Data/Analysis/8_OTUs_beta/$dir2 2  ;	
done

for i in $(ls $Data/14_Statistics/beta/*/*.otu_table*.txt);
do
	dir1=${i#*beta}
	dir2=`echo "$dir1" | cut -d / -f 2`
	[ ! -d $Data/Analysis/8_OTUs_beta/$dir2 ] && mkdir -p $Data/Analysis/8_OTUs_beta/$dir2;
#	cp $i $Data/Analysis/8_OTUs_beta/$dir2;
	Rscript $Data/bin/A7.2_Plot_OTUs_beta_bar.r $i $Data/bin/$dir2-color.txt $Data/Analysis/8_OTUs_beta/$dir2 2  ;	

done
