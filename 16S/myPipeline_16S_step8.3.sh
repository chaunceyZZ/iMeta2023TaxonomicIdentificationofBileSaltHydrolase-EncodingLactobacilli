#StepIV:Alignment and remove Chimera
#NOTE: The QIIME enviroment should be set firstly!
#### docker run -it -v `pwd`:/home --name=xxxx fischuu/qiime-1.9.1
Data=$(dirname $PWD)
# GREENGENE=/home/Project/Ref_database/greengenes;
GREENGENE=$Data/bin/greengenes


####
###	sub3W1 :note modify
####


###
###	gut-grouping:#SampleID Group
###


###
###	gr=(Total,MC,MeS,MeH,MA,MAS) :note modify
###

############	Note:Less than SubsampleSeq,if too large,will be missing files in 12_alphaRarefactionPlot;less 1W
##### Note
###	Get min
#	for i in $(ls $Data/3_filtered_fa/*.fa);do cat $i | grep -E ">" | wc -l;done

###	for i in $(ls $Data/3_filtered_fa/*.fa);do cat $i | grep -E ">" | wc -l;done

# SubsampleSeq=`(for i in $(ls $Data/3_filtered_fa/*.fna);do cat $i | grep -E ">" | wc -l;done) | sort  -n | head -1`
# echo $SubsampleSeq
# 
# Depth_alpha=$SubsampleSeq
# Depth_beta=$SubsampleSeq

Depth_alpha=40000
Depth_beta=40000
############

########### NOTE
gr=("Total")
###########



## Beta diversity
echo "13_betaDiv"
[ ! -d $Data/13_betaDiv ] && mkdir $Data/13_betaDiv

echo "single_rarefaction"
for i in $(ls $Data/7_Taxa_Greengene/*.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom)
do
	ID1=$(basename ${i%.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom})
	ID2=$(basename ${i%.biom})_even${Depth_beta}.biom
	[ ! -d $Data/13_betaDiv/$ID1 ] && mkdir $Data/13_betaDiv/$ID1
	single_rarefaction.py -i $i -o $Data/13_betaDiv/$ID1/$ID2 -d $Depth_beta;
done

echo "beta_diversity"
for i in $(ls $Data/13_betaDiv/*/*.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table_even${Depth_beta}.biom)
do
	ID1=$(basename ${i%.gut*})
	ID2=$(basename ${i%.filter*})_rep_set.filter_4_noChimera_tree.tre
	beta_diversity.py -i $i -o $Data/13_betaDiv/$ID1/ -t $Data/6_Phylogeny_Greengene/$ID2 -m unifrac,unifrac_g,unweighted_unifrac,weighted_unifrac;
done

echo "principal_coordinates"
for m in unifrac unweighted_unifrac weighted_unifrac;
do 
	for i in ${gr[@]}
	do
		principal_coordinates.py -i $Data/13_betaDiv/$i/${m}_$i.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table_even${Depth_beta}.txt -o $Data/13_betaDiv/$i/${m}_$i.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table_even${Depth_beta}.pc.txt;
	done
done

echo "complete!!!"
