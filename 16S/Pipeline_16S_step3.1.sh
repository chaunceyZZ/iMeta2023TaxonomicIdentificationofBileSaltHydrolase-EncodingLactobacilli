#StepII:Filter the raw reads by the minimum length 400bp and subsample 40000 reads from each sample!
#NOTE: The QIIME enviroment should be set firstly!
### 安装 400 和 subsample过滤数据

Data=$(dirname $PWD)
cd $Data;
[ ! -d $Data/3_filtered_fa ] && mkdir $Data/3_filtered_fa;

#2.1 filter Fasta by minimum length 400bp
SIZE=400;  # v3-v5
###	Note: *.fa
###	Change
#	for i in $(ls $Data/0_raw_data/*.fa);
#	do
	
#	mv $i $(dirname $i)/$(basename ${i%.fna}).fa
#	done


for fa in $(ls $Data/2_combine_fa/*.fna);
do
	$Data/bin/filterFasta_byLength.pl $fa $SIZE; ##### 以size过滤
done
mv  $Data/2_combine_fa/*.filter_$SIZE.fna $Data/3_filtered_fa;

#2.2 subsample 40000 reads from each sample!
##### Note
###	Get min
#	for i in $(ls $Data/2_combine_fa/*.fa);
#	do
#		cat $i | grep -E ">" | wc -l
#	done

#  for i in $(ls *.fna); do cat $i | grep -E ">" | wc -l; done

