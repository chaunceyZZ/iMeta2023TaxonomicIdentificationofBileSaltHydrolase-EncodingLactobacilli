#Join the pair_end reads to one longer reads!
#Usually, the Sequencing company will do this step!

### 双端序列合并
Data=$(dirname $PWD);

[ ! -d $Data/1_combine_fq ] && mkdir $Data/1_combine_fq;
## NOTE: *.1.fq
cd $Data/1_combine_fq
:>commandlist_join_paired_ends
for i in $(ls $Data/0_raw_data/*.1.fq);
do
	echo $(basename $i);
	echo "join_paired_ends.py -f $i -r ${i%.1.fq}.2.fq -o ${i%.1.fq}.join;mv ${i%.1.fq}.join/fastqjoin.join.fastq $Data/1_combine_fq/$(basename $i .1.fq).join.fq;rm -r ${i%.1.fq}.join;" >> commandlist_join_paired_ends
done

ParaFly -c commandlist_join_paired_ends -CPU 30


##join_paired_ends.py
##双端序列合并
##任务是把双端序列合并，根据两端序列末端的互补配对，可以合变为我们扩增区域的序列，同时还可以对重叠区的质量进行校正，保留最高测序质量的碱基结果


