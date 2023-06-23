#Analysis StepVIII: do PICRUSt to predict functional profiling
#NOTE: The QIIME enviroment should be set firstly!

Data=$(dirname $PWD)
# GREENGENE=/home/Project/Ref_database/greengenes;
GREENGENE=$Data/bin/greengenes

cd $Data;
[ ! -d $Data/Analysis/9_PICRUSt ] && mkdir $Data/Analysis/9_PICRUSt;

###
###	NOTE:sub3W
###

##	CP   Total.gut.sub3W.fa to $Data/Analysis/9_PICRUSt

for i in $(ls $Data/3_filtered_fa/*.gut.sub3W1.fa)
do
	ff=$(basename $i)
	f1=`echo "$ff" | cut -d . -f 1`
	[ ! -d $Data/Analysis/9_PICRUSt/$f1 ] && mkdir $Data/Analysis/9_PICRUSt/$f1;
	cp $i $Data/Analysis/9_PICRUSt/$f1;

done


#Step 8.1 Pick OTUs
for i in $(ls $Data/Analysis/9_PICRUSt/*/*.gut.sub3W1.fa)
do
	dir=$(dirname $i)
	ff=$(basename $i)
	f1=`echo "$ff" | cut -d . -f 1`
	f2=$f1".otu_table.biom"
	pick_closed_reference_otus.py -f -i $i -o $dir -p $Data/bin/otu_picking_params_97.txt -r $GREENGENE/gg_13_5_otus/rep_set/97_otus.fasta -t $GREENGENE/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt
	mv $dir/otu_table.biom $dir/$f2

done



#Step 8.2 Normalize OTUs
for i in $(ls $Data/Analysis/9_PICRUSt/*/*.otu_table.biom)
do
	dir=$(dirname $i)
	ff=$(basename $i)
	f1=`echo "$ff" | cut -d . -f 1`
	f2=$f1".normolized.out.biom"

	
	normalize_by_copy_number.py -i $i -o $dir/$f2

done



#Step 8.3 Predict KEGG
for i in $(ls $Data/Analysis/9_PICRUSt/*/*.normolized.out.biom)
do

	dir=$(dirname $i)
	ff=$(basename $i)
	f1=`echo "$ff" | cut -d . -f 1`

	f2=$f1".predictions.kegg.tab"
	f3=$f1".predictions.kegg.biom"
	f4=$f1".nsti_per_sample.kegg.tab"

	predict_metagenomes.py -f -i $i -o $dir/$f2
	predict_metagenomes.py -i $i -o $dir/$f3 -a $dir/$f4

done






#Step 8.4 Predict COG
for i in $(ls $Data/Analysis/9_PICRUSt/*/*.normolized.out.biom)
do
	dir=$(dirname $i)
	ff=$(basename $i)
	f1=`echo "$ff" | cut -d . -f 1`

	f2=$f1".predictions.cog.tab"
	f3=$f1".predictions.cog.biom"
	f4=$f1".nsti_per_sample.cog.tab"

	predict_metagenomes.py --type_of_prediction cog -f -i $i -o $dir/$f2
	predict_metagenomes.py --type_of_prediction cog -i $i -o $dir/$f3 -a $dir/$f4

done


#Step 8.5 Predict rfam
for i in $(ls $Data/Analysis/9_PICRUSt/*/*.normolized.out.biom)
do
	dir=$(dirname $i)
	ff=$(basename $i)
	f1=`echo "$ff" | cut -d . -f 1`

	f2=$f1".predictions.rfam.tab"
	f3=$f1".predictions.rfam.biom"
	f4=$f1".nsti_per_sample.rfam.tab"

	predict_metagenomes.py --type_of_prediction rfam -f -i $i -o $dir/$f2
	predict_metagenomes.py --type_of_prediction rfam -i $i -o $dir/$f3 -a $dir/$f4

done




#Step 8.6 Collapse KEGG table data 
for i in $(ls $Data/Analysis/9_PICRUSt/*/*.predictions.kegg.biom)
do
	dir=$(dirname $i)
	ff=$(basename $i)
	f1=`echo "$ff" | cut -d . -f 1`

	f2=$f1".catagrized.kegg.L2.biom"
	f3=$f1".catagrized.kegg.L2.tab"
	
	categorize_by_function.py -i $i -c KEGG_Pathways -l 2 -o $dir/$f2
	categorize_by_function.py -f -i $i -c KEGG_Pathways -l 2 -o $dir/$f3

done




#Step 8.7 Collapse COG table data
for i in $(ls $Data/Analysis/9_PICRUSt/*/*.predictions.cog.biom)
do
	dir=$(dirname $i)
	ff=$(basename $i)
	f1=`echo "$ff" | cut -d . -f 1`

	f2=$f1".catagrized.cog.L2.biom"
	f3=$f1".catagrized.cog.L2.tab"

	categorize_by_function.py -i $i -c COG_Category -l 2 -o $dir/$f2
	categorize_by_function.py -f -i $i -c COG_Category -l 2 -o $dir/$f3


done




#Step 8.8 Summrize Plot
for i in $(ls $Data/Analysis/9_PICRUSt/*/*.catagrized.kegg.L2.biom)
do
	dir=$(dirname $i)
	ff=$(basename $i)
	f1=`echo "$ff" | cut -d . -f 1`

	f2="KEGG_"$f1

	summarize_taxa_through_plots.py -i $i  -p $Data/bin/qiime.kegg.params.txt -o $dir/$f2
done




for i in $(ls $Data/Analysis/9_PICRUSt/*/*.catagrized.cog.L2.biom)
do
	dir=$(dirname $i)
	ff=$(basename $i)
	f1=`echo "$ff" | cut -d . -f 1`

	f2="COG_"$f1

	summarize_taxa_through_plots.py -i $i  -p $Data/bin/qiime.cog.params.txt -o $dir/$f2

done

#Step 8.9 Count to Abundance
for i in $(ls $Data/Analysis/9_PICRUSt/*/*L2.tab);do sed -i '1d' $i;done
for i in $(ls $Data/Analysis/9_PICRUSt/*/*L2.tab);do perl -pi -e "s/\#//g" $i;done
#####################

for i in $(ls $Data/Analysis/9_PICRUSt/*/*L2.tab);
do
	dd=$(dirname $i)

	ff=$(basename $i)
	g1=`echo "$ff" | cut -d . -f 1`
	
	
	Rscript $Data/bin/A8.1_PICRUSt_Count2Abundance.r $i $Data/bin/$g1-color.txt $dd 2;

done
for i in $(ls $Data/Analysis/9_PICRUSt/*/*abundance.txt);
do
	perl -pi -e "s/^\t/ID\t/g" $i;
done

#############################

