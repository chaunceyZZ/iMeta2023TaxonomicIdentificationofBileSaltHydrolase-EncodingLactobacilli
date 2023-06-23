#StepIV:Alignment and remove Chimera
#NOTE: The QIIME enviroment should be set firstly!
######### 注释库为 greenenes 
Data=$(dirname $PWD)
# GREENGENE=/home/Project/Ref_database/greengenes;
GREENGENE=$Data/bin/greengenes

[ ! -d $Data/6_Phylogeny_Greengene ] && mkdir $Data/6_Phylogeny_Greengene

####
###	sub3W1 :note modify
####


#Step 6.1 Make Phylogeny
echo "6.1"
for i in $(ls $Data/5_remove_Chimera_Greengene/*.gut.sub3W1_rep_set.filter_4_noChimera_aligned_pfiltered.fasta)
do
	ID1=$(basename ${i%_aligned_pfiltered.fasta})_tree.tre
	make_phylogeny.py -i $i -o $Data/6_Phylogeny_Greengene/$ID1;
	###通过把代表otus的序列进行多序列比对构建进化树
done
#Step 6.2 Taxanomy Assign
echo "6.2"
[ ! -d $Data/7_Taxa_Greengene ] && mkdir $Data/7_Taxa_Greengene;
for i in $(ls $Data/4_Picked_OTUs/*.gut.sub3W1_rep_set.filter_4.fa)
do
	assign_taxonomy.py -i $i -o $Data/7_Taxa_Greengene -t $GREENGENE/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt -r $GREENGENE/gg_13_8_otus/rep_set/97_otus.fasta;
	########将分类标准分配给OTU代表序列
done
#Step 6.3 Make OTUs table
echo "6.3"
for i in $(ls $Data/4_Picked_OTUs/*.gut.sub3W1_otus.filter_4.txt)
do
	ID1=$(basename ${i%_otus.filter_4.txt})_rep_set.filter_4_tax_assignments.txt
	ID2=$(basename ${i%_otus.filter_4.txt}).filter_4_noChimera_rmAlignFail.otu_table.biom
	ID3=$(basename ${i%_otus.filter_4.txt}).rmOTUs_Chimera_AlignFail.txt
	make_otu_table.py -i $i -t $Data/7_Taxa_Greengene/$ID1 -o $Data/7_Taxa_Greengene/$ID2 -e $Data/5_remove_Chimera_Greengene/$ID3;
	#####  制作OTU表
done
#Step 6.4 Convert OTUs table from biom-format to txt-format
echo "6.4"
for i in $(ls $Data/7_Taxa_Greengene/*.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom)
do
	ID1=$(basename ${i%.biom}).txt
	biom convert -i $i -o $Data/7_Taxa_Greengene/$ID1 --table-type "OTU table" --to-tsv --header-key="taxonomy" --output-metadata-id="Consensus Lineage"
done
#Step 6.5 Summarize Taxanomy
echo "6.5"
for i in $(ls $Data/7_Taxa_Greengene/*.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom)
do
	summarize_taxa.py -i $i  -o $Data/7_Taxa_Greengene/
	ID1=$(basename ${i%.biom}).summary.txt
	biom summarize-table -i $i -o $Data/7_Taxa_Greengene/$ID1

done


for L in {2..6};
do 
	for i in $(ls $Data/7_Taxa_Greengene/*.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table_L${L}.biom)
	do
		ID1=$(basename ${i%.biom}).summary.txt
		biom summarize-table -i $i -o $Data/7_Taxa_Greengene/$ID1;
	done

done


