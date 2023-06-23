Data=$(dirname $PWD)
cd $Data;
SIZE=400;
SubsampleSeq=56000; # 56248
##############	NOTE:sub3W1

# SubsampleSeq=`(for i in $(ls $Data/3_filtered_fa/*.fna);do cat $i | grep -E ">" | wc -l;done) | sort  -n | head -1`
echo $SubsampleSeq

###
for fa in $(ls $Data/3_filtered_fa/*.filter_$SIZE.fna)
do
	seq_num=$(grep -c ">" $fa)
	perc=$(echo "$SubsampleSeq $seq_num"| awk -F' ' '{printf "%.3f",$1/$2}')
	subsample_fasta.py -i $fa -p $perc -o ${fa%.filter*}.sub3W1.fa    ####随机抽取指定数目的FASTA序列

done
#############	NOTE:Modify the original file	################
# 参数1为输出文件名;2，3为输入文件;4为参数2,3是否一致，1表示一致，2表示不一致;5为参数,2表示Pair,1表示Single
###	Paired
perl $Data/bin/merge-indi-fasta.pl $Data/3_filtered_fa/Total.gut.sub3W1.fa "$Data/3_filtered_fa/*.join.sub3W1.fa" "$Data/3_filtered_fa/*.join.sub3W1.fa" 1 2

#### Single
#perl $Data/bin/merge-indi-fasta.pl $Data/3_filtered_fa/Total.gut.sub3W1.fa "$Data/3_filtered_fa/*.sub3W1.fa" "$Data/3_filtered_fa/*.sub3W1.fa" 1 1



#perl $Data/bin/merge-indi-fasta.pl $Data/3_filtered_fa/MC.gut.sub3W1.fa "$Data/3_filtered_fa/Mdl*.join.sub3W1.fa" "$Data/3_filtered_fa/Ctrl*.join.sub3W1.fa" 2 2
#perl $Data/bin/merge-indi-fasta.pl $Data/3_filtered_fa/MeS.gut.sub3W1.fa "$Data/3_filtered_fa/Met*.join.sub3W1.fa" "$Data/3_filtered_fa/SERP*.join.sub3W1.fa" 2 2
#perl $Data/bin/merge-indi-fasta.pl $Data/3_filtered_fa/MeH.gut.sub3W1.fa "$Data/3_filtered_fa/Met*.join.sub3W1.fa" "$Data/3_filtered_fa/HQ*.join.sub3W1.fa" 2 2
#perl $Data/bin/merge-indi-fasta.pl $Data/3_filtered_fa/MA.gut.sub3W1.fa "$Data/3_filtered_fa/Mdl[0-9]*.join.sub3W1.fa" "$Data/3_filtered_fa/Ab[0-9]*.join.sub3W1.fa" 2
#perl $Data/bin/merge-indi-fasta.pl $Data/3_filtered_fa/MAS.gut.sub3W1.fa "$Data/3_filtered_fa/Mdl*.join.sub3W1.fa" "$Data/3_filtered_fa/AbS*.join.sub3W1.fa" 2 2

