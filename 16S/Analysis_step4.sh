#Analysis StepI: Classify the abundance data 
#Current working taxanomic level is Genus

Data=$(dirname $PWD)

[ ! -d $Data/Analysis ] && mkdir $Data/Analysis;
[ ! -d $Data/Analysis/4_common_OTUs ] && mkdir $Data/Analysis/4_common_OTUs;

for i in $(ls $Data/7_Taxa_Greengene/*.otu_table.txt)
do
	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	[ ! -d $Data/Analysis/4_common_OTUs/$c1 ] && mkdir $Data/Analysis/4_common_OTUs/$c1;
	cp $i $Data/Analysis/4_common_OTUs/$c1/;
done


for i in $(ls $Data/Analysis/4_common_OTUs/*/*.otu_table.txt)
do
	perl $Data/bin/A4.1_extract_OTUs_taxa.pl $i;
done

for i in $(ls $Data/Analysis/4_common_OTUs/*/*.otu.counts.txt)
do
	perl $Data/bin/A4.2_calculate_OTUs_abundance.pl $i;
done


[ ! -d $Data/Analysis/4_common_OTUs/Total/Total ] && mkdir $Data/Analysis/4_common_OTUs/Total/Total

for i in $(ls $Data/Analysis/4_common_OTUs/*/*.otu.*.txt)
do
	perl $Data/bin/A4.3_extract_commonOTUs.pl $i;
done

for i in $(ls $Data/Analysis/4_common_OTUs/*/*.otu.*.txt)
do
	ID=$(basename $i)
	echo $ID
	c1=`echo "$ID" | cut -d . -f 1`
	perl $Data/bin/A4.4_extract_commonOTUs-groups.pl $i $Data/bin/$c1-color.txt;
done

for i in $(ls $Data/Analysis/4_common_OTUs/*/*.otu.counts.txt)
do
	ID=$(basename $i)
	c1=`echo "$ID" | cut -d . -f 1`
	Rscript $Data/bin/A4.5_Plot_CommonOTUs.r $i $Data/Analysis/4_common_OTUs/$c1/;
done

