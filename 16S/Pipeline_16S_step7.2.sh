#StepIV:Alignment and remove Chimera
#NOTE: The QIIME enviroment should be set firstly!

Data=$(dirname $PWD)
# GREENGENE=/home/Project/Ref_database/greengenes;
GREENGENE=$Data/bin/greengenes


####
###	sub3W1 :note modify
####

### gut-grouping:#SampleID Group

#Step 7.2: Make OTU network, visulizing by Cytoscape
echo "7.2"
for i in $(ls $Data/7_Taxa_Greengene/*.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom)
do
	ID1=$(basename ${i%.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom})-grouping.info
	make_otu_network.py -i $i -m $Data/bin/$ID1 -o $Data/8_Plot -b Group,"Group";
done
###
###	gut-grouping.info :#
###



#Step 7.3: Make OTU heatmap
echo "7.3"
for i in $(ls $Data/7_Taxa_Greengene/*.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom)
do
	ID1=$(basename ${i%.biom})_heatmap.pdf
	
	make_otu_heatmap.py -i $i -o $Data/8_Plot/$ID1  --width 7 --height 7
	ID2=$(basename ${i%.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom})-grouping.info
	make_otu_heatmap.py -i $i -o $Data/8_Plot/$ID1 -c "Group" -m $Data/bin/$ID2 --width 7 --height 7
done


