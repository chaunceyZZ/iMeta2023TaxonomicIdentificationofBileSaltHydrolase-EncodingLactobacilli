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
SubsampleSeq=56000
echo $SubsampleSeq

Depth_alpha=$SubsampleSeq
Depth_beta=$SubsampleSeq
# Depth_alpha=23000
# Depth_beta=23000
############

########### NOTE
gr=("Total")
###########


#Step 8.1: Multiple rarefaction and alpha diversity
echo "8.1"
[ ! -d $Data/9_Rarefaction ] && mkdir $Data/9_Rarefaction
[ ! -d $Data/10_alphaDiv ] && mkdir $Data/10_alphaDiv

for i in $(ls $Data/7_Taxa_Greengene/*.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom)
do
	ID1=$(basename ${i%.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table.biom})
	multiple_rarefactions.py -i $i -m 100 -x $Depth_alpha -s 100 -n 10 -o $Data/9_Rarefaction/$ID1/rare100_$Depth_alpha/
	echo "multiple_rarefactions.py complete!" 
	ID2=$(basename ${i%.filter_4_noChimera_rmAlignFail.otu_table.biom})_rep_set.filter_4_noChimera_tree.tre
	alpha_diversity.py -i $Data/9_Rarefaction/$ID1/rare100_$Depth_alpha/ -m shannon,PD_whole_tree,chao1,observed_otus -o $Data/10_alphaDiv/$ID1 -t $Data/6_Phylogeny_Greengene/$ID2
	echo "alpha_diversity.py com!"

	for j in {0..9};do mkdir $Data/9_Rarefaction/$ID1/rare100_$Depth_alpha/iter_${j}; mv $Data/9_Rarefaction/$ID1/rare100_$Depth_alpha/rarefaction_*_${j}.biom  $Data/9_Rarefaction/$ID1/rare100_$Depth_alpha/iter_${j}/ ;done

	for m in {0..9};do mkdir $Data/10_alphaDiv/$ID1/iter_${m};done

	for n in {0..9};do alpha_diversity.py -i $Data/9_Rarefaction/$ID1/rare100_$Depth_alpha/iter_${n}/ -m shannon,PD_whole_tree,chao1,observed_otus -o $Data/10_alphaDiv/$ID1/iter_${n}/ -t $Data/6_Phylogeny_Greengene/$ID2;done

	for k in {0..9};do mv $Data/10_alphaDiv/$ID1/iter_${k}/* $Data/10_alphaDiv/$ID1; rmdir $Data/10_alphaDiv/$ID1/iter_${k};done

done



echo "11_collated_alpha"
[ ! -d $Data/11_collated_alpha ] && mkdir $Data/11_collated_alpha


for L in ${gr[@]}
do
	[ ! -d $Data/11_collated_alpha/$L ] && mkdir $Data/11_collated_alpha/$L
	collate_alpha.py -i $Data/10_alphaDiv/$L/ -o $Data/11_collated_alpha/$L/
done
