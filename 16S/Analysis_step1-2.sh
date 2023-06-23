#!/usr/bin/env bash

database=/home/zhouxingchen/Downloads/database_16S/minikraken_8GB_20200312

Data=$(dirname $PWD)

[ ! -d $Data/Analysis ] && mkdir $Data/Analysis;
[ ! -d $Data/Analysis/0_absolute_Abundance ] && mkdir -p $Data/Analysis/0_absolute_Abundance;

# :>$Data/Analysis/0_absolute_Abundance/Total.asvclass.txt
# for i in $(ls $Data/2_combine_fa/*.join.fna);
# do
# 	kraken2 --db $database \
# 		--threads 80 \
# 		--report $Data/Analysis/0_absolute_Abundance/$(basename $i .join.fna).kreport2  \
# 		--output $Data/Analysis/0_absolute_Abundance/$(basename $i .join.fna).result2 \
# 		--unclassified-out $Data/Analysis/0_absolute_Abundance/$(basename $i .join.fna).unclassified.fa \
# 		--classified-out $Data/Analysis/0_absolute_Abundance/$(basename $i .join.fna).classified.fa  \
# 		$i  2>&1 | \
# 		tee $Data/Analysis/0_absolute_Abundance/$(basename $i .join.fna).log
# 		# --use-mpa-style \
# 		# $i
# 
# 	awk '$1 == "C" && $3 != 0' $Data/Analysis/0_absolute_Abundance/$(basename $i .join.fna).result2 \
# 		| awk '{$0=$2 "\t"$3}1' \
# 		| sed '1i\asvid\tspid' > $Data/Analysis/0_absolute_Abundance/$(basename $i .join.fna).asvclass.txt
# 
# 	python mappingArea.py $i $Data/Analysis/0_absolute_Abundance/$(basename $i .join.fna).csv
# 
#         awk '$1 == "C" && $3 != 0' $Data/Analysis/0_absolute_Abundance/$(basename $i .join.fna).result2 \
#                 | awk -v A=$(basename $i .join.fna) '{$0=A"\t"$3}1' >> $Data/Analysis/0_absolute_Abundance/Total.asvclass.txt
# done
# 
# sed -i '1i\sample\tspid' $Data/Analysis/0_absolute_Abundance/Total.asvclass.txt



########################### pre abundance ####################################
# /home/zhouxingchen/miniconda3/bin/Rscript A0.1_absolute_abundance_pre.R $Data/Analysis/0_absolute_Abundance /home/zhouxingchen/Downloads/db/taxonomy/merged.dmp $Data/Analysis/0_absolute_Abundance/Total.asvclass.txt /home/zhouxingchen/Downloads/db/kraken2/dbpy.txt
# 
# 
# 
# for i in $(ls $Data/Analysis/0_absolute_Abundance2/*.PlyGen.absolute.abundance.txt);
# do
#         ID=$(basename $i);
#         Rscript $Data/bin/A1.1_Plot_Abundance_heatmap.r $i $Data/Analysis/${ID%.Ply*}-color3.txt $Data/Analysis/0_absolute_Abundance2 $Data/bin;
# done

# for i in $(ls $Data/Analysis/0_absolute_Abundance2/*.Genus.absolute.abundance.txt);
# do
#         ID=$(basename $i);
#         Rscript $Data/bin/A1.2_Plot_Abundance_box.r $i $Data/Analysis/${ID%.Genus*}-color3.txt $Data/Analysis/0_absolute_Abundance2 15
#         Rscript $Data/bin/A1.2.2_Plot_Abundance_between_box.r $i $Data/Analysis/${ID%.Genus*}-color3.txt $Data/Analysis/0_absolute_Abundance2 15 "L.animalis-Control" "L.fermentum-Control" "L.salvarius-Control" "L.reuteri-Control"
# done;

for i in $(ls $Data/Analysis/0_absolute_Abundance/*.PlyGen.absolute.abundance.txt);
do
        Out_dir=$(dirname $i)
        File=$(basename $i)
        p=${File/PlyGen/Phylum};
        echo $p
        c1=`echo "$File" | cut -d . -f 1`
        Rscript $Data/bin/A1.4_Abundance_bar.r $Out_dir $i $p $Data/Analysis/$c1-color4.txt

done
