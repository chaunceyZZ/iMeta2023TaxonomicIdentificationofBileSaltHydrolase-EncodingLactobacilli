#Analysis StepI:Convert the raw Qiime abundance file to the tab-format files!
#Current working taxanomic level is Genus

Data=$(dirname $PWD)

[ ! -d $Data/Analysis ] && mkdir $Data/Analysis;
[ ! -d $Data/Analysis/1_Abundance2 ] && mkdir $Data/Analysis/1_Abundance2;

# cp $Data/7_Taxa_Greengene/*.otu_table_L6.txt $Data/Analysis/1_Abundance;
# for i in $(ls $Data/Analysis/1_Abundance/*.otu_table_L6.txt);
# do
# 	perl $Data/bin/A0.0_extract_abundance.pl $i;
# done
# 
# for i in $(ls $Data/Analysis/1_Abundance/*.abundance.tmp1.txt);
# do
# 	perl $Data/bin/A0.1_rm_duplicate_taxa.pl $i;
# done
# 
# for i in $(ls $Data/Analysis/1_Abundance/*.abundance.tmp2.txt);
# do
# 	perl $Data/bin/A0.2_sort_indi_abundance.pl $i;
# done
# 
# rm $Data/Analysis/1_Abundance/*abundance.tmp*.txt;
# 
# for i in $(ls $Data/Analysis/1_Abundance/*.otu_table_L6.txt);
# do
# 	perl $Data/bin/A0.3_add_phyla_abundance.pl $i;
# done
# for i in $(ls $Data/Analysis/1_Abundance/*.otu_table_L6.txt);
# do
# 	perl $Data/bin/A0.0_extract_abundance_phylum.pl $i;
# done


# for i in $(ls $Data/Analysis/1_Abundance2/*.PlyGen.abundance.txt);
# do
# 	ID=$(basename $i);
# 	Rscript $Data/bin/A1.1_Plot_Abundance_heatmap.r $i $Data/Analysis/${ID%.Ply*}-color2.txt $Data/Analysis/1_Abundance2 $Data/bin;
# done
# 
# for i in $(ls $Data/Analysis/1_Abundance2/*.Genus.abundance.txt);
# do
# 	ID=$(basename $i);
# 	Rscript $Data/bin/A1.2_Plot_Abundance_box.r $i $Data/Analysis/${ID%.Genus*}-color2.txt $Data/Analysis/1_Abundance2 15
# 	Rscript $Data/bin/A1.2.2_Plot_Abundance_between_box.r $i $Data/Analysis/${ID%.Genus*}-color2.txt $Data/Analysis/1_Abundance2 15 "L.animalis-Control" "L.fermentum-Control" "L.salvarius-Control" "L.reuteri-Control" 
# done;


for i in $(ls $Data/Analysis/1_Abundance0412/*.PlyGen.abundance.txt);
do
	Out_dir=$(dirname $i)
	File=$(basename $i)
	p=${File/PlyGen/Phylum};
	echo $p
	c1=`echo "$File" | cut -d . -f 1`
	Rscript $Data/bin/A1.4_Abundance_bar.r $Out_dir $i $p $Data/Analysis/$c1-color3.txt

done


###Group###
#bash A1.0_Extract_Group_abundance.sh
#bash A1.3_Group_figures.sh

