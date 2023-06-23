#StepVIII:Doing some statistics
#NOTE: The QIIME enviroment should be set firstly!

Data=$(dirname $PWD)

####
###	sub3W1 :note modify
####

############	Note:Less than SubsampleSeq
Depth_alpha=40000
Depth_beta=40000
############


[ ! -d $Data/14_Statistics ] && mkdir $Data/14_Statistics;

#####	Alpha-diversity
echo "Alpha-diversity"
for i in $(ls $Data/7_Taxa_Greengene/*.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom)
do
	ID0=$(basename ${i%.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom})
	[ ! -d $Data/14_Statistics/alpha/$ID0 ] && mkdir -p $Data/14_Statistics/alpha/$ID0

	echo "Step 14.1 Calculate GOOD COVERAGE (Alpha-diversity)!"
	ID1=$(basename ${i%.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom}).final.otu.good.coverage.txt
	alpha_diversity.py -i $i -m goods_coverage -o $Data/14_Statistics/alpha/$ID0/$ID1;

	echo "Step 14.2 Calculate Shanno Index (Alpha-diversity)!"
	ID2=$(basename ${i%.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom}).final.otu.shannon.txt
	alpha_diversity.py -i $i -m shannon -o $Data/14_Statistics/alpha/$ID0/$ID2;

	echo "Step 14.3 Calculate Chao1 (Alpha-diversity)!"
	ID3=$(basename ${i%.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom}).final.otu.chao1.txt
	alpha_diversity.py -i $i -m chao1 -o $Data/14_Statistics/alpha/$ID0/$ID3;

	echo "Step 14.4 Calculate observed_species (Alpha-diversity)!"
	ID4=$(basename ${i%.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom}).final.otu.observed_species.txt
	alpha_diversity.py -i $i -m observed_species -o $Data/14_Statistics/alpha/$ID0/$ID4;

	echo "Step 14.5 Calculate simpson Index (Alpha-diversity)!"
	ID5=$(basename ${i%.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom}).final.otu.simpson.txt
	alpha_diversity.py -i $i -m simpson -o $Data/14_Statistics/alpha/$ID0/$ID5;
done

####	Beta-diversity
echo "Beta-diversity"
for i in $(ls $Data/13_betaDiv/*/*.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table_even${Depth_beta}.biom)
do
	ID0=$(basename ${i%.gut*})
	[ ! -d $Data/14_Statistics/beta/$ID0 ] && mkdir -p $Data/14_Statistics/beta/$ID0

	echo "Step 14.6 Calculate bray_curtis Index  (Beta-diversity)!"
	beta_diversity.py -i $i -m bray_curtis -o $Data/14_Statistics/beta/$ID0/

	echo "Step 14.7 Calculate bray_curtis_magurran Index  (Beta-diversity)!"
	beta_diversity.py -i $i -m bray_curtis_magurran -o $Data/14_Statistics/beta/$ID0/

	echo "Step 14.8 Calculate binary_jaccard Index  (Beta-diversity)!"
	beta_diversity.py -i $i -m binary_jaccard -o $Data/14_Statistics/beta/$ID0/

	echo "Step 14.9 Calculate abund_jaccard Index  (Beta-diversity)!"
	beta_diversity.py -i $i -m abund_jaccard -o $Data/14_Statistics/beta/$ID0/

	echo "Step 14.10 Calculate spearman_approx Index  (Beta-diversity)!"
	beta_diversity.py -i $i -m spearman_approx -o $Data/14_Statistics/beta/$ID0/

	echo "Step 14.11 Calculate euclidean Index  (Beta-diversity)!"
	beta_diversity.py -i $i -m euclidean -o $Data/14_Statistics/beta/$ID0/

	
done
  
echo "complete!"
