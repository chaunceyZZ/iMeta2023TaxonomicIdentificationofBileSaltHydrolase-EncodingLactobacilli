#StepIV:Alignment and remove Chimera
#NOTE: The QIIME enviroment should be set firstly!

Data=$(dirname $PWD)
# GREENGENE=/home/Project/Ref_database/greengenes;
GREENGENE=$Data/bin/greengenes


####
###	sub3W1 :note modify
####

### gut-grouping:#SampleID Group

#Step 7.1: Plot Taxa distributions
echo "7.1"
[ ! -d $Data/8_Plot ] && mkdir $Data/8_Plot
for i in $(ls $Data/7_Taxa_Greengene/*.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table_L2.txt)
do
	ID1=$(basename ${i%_L2.txt})_L3.txt
	ID2=$(basename ${i%_L2.txt})_L4.txt
	ID3=$(basename ${i%_L2.txt})_L5.txt
	ID4=$(basename ${i%_L2.txt})_L6.txt
	plot_taxa_summary.py -i $i,$Data/7_Taxa_Greengene/$ID1,$Data/7_Taxa_Greengene/$ID2,$Data/7_Taxa_Greengene/$ID3,$Data/7_Taxa_Greengene/$ID4 -l Phylum,Class,Order,Family,Genus -o $Data/8_Plot;
done

