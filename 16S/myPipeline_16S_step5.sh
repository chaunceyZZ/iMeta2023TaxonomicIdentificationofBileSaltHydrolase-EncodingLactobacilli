#StepIV:Alignment and remove Chimera
#NOTE: The QIIME enviroment should be set firstly!
##### 对齐删除冗余
Data=$(dirname $PWD)
# GREENGENE=/home/Project/Ref_database/greengenes;
GREENGENE=$Data/bin/greengenes

[ ! -d $Data/5_remove_Chimera_Greengene ] && mkdir $Data/5_remove_Chimera_Greengene;

####
###	sub3W1 :note modify
####


#Step 4.1
echo "4.1"
for i in $(ls $Data/4_Picked_OTUs/*gut.sub3W1_rep_set.filter_4.fa)
do
	align_seqs.py -i $i -m pynast -o $Data/5_remove_Chimera_Greengene -t $GREENGENE/core_set_aligned.fasta.imputed
	####根据选择的方法来输入序列之间的比对以及与参考序列的比对
done
#Step 4.2
echo "4.2"
for i in $(ls $Data/5_remove_Chimera_Greengene/*.gut.sub3W1_rep_set.filter_4_aligned.fasta)
do
	ID1=$(basename ${i%.fasta}).chimeric_seqs.txt
	identify_chimeric_seqs.py -i $i -m ChimeraSlayer -a $GREENGENE/core_set_aligned.fasta.imputed -o $Data/5_remove_Chimera_Greengene/$ID1;
	####去除嵌合体（PCR扩增过程中产生的序列是来自多个模版）
done
#Step 4.3
echo "4.3"
for i in $(ls $Data/5_remove_Chimera_Greengene/*.gut.sub3W1_rep_set.filter_4_aligned.fasta)
do
	ID1=$(basename ${i%_aligned.fasta})_noChimera_aligned.fasta
	ID2=$(basename ${i%.fasta}).chimeric_seqs.txt
	filter_fasta.py -f $i -o $Data/5_remove_Chimera_Greengene/$ID1 -s $Data/5_remove_Chimera_Greengene/$ID2 -n;
	###从fasta中提取或者过滤掉多个序列
done
#Step 4.4
echo "4.4"
for i in $(ls $Data/5_remove_Chimera_Greengene/*.gut.sub3W1_rep_set.filter_4_noChimera_aligned.fasta)
do
	filter_alignment.py -i $i -o $Data/5_remove_Chimera_Greengene;
	### 过滤比对后的数据
	###有两点作用：1它是去掉在所有序列中都是gap的位置（特别是PyNAST，200-400个碱基的序列，跟16S全长基因比对后的序列）；2，同时可以提供lanemask，指定那些位置在构建进化树的时候需要保留，哪一些需要去掉，这是针对保守位点的需要。
done
#Step 4.5
echo "4.5"
for i in $(ls $Data/5_remove_Chimera_Greengene/*.gut.sub3W1_rep_set.filter_4_failures.fasta)
do
	ID1=$(basename ${i%_failures.fasta})_aligned.chimeric_seqs.txt
	ID2=$(basename ${i%_rep_set.filter_4_failures.fasta}).rmOTUs_Chimera_AlignFail.txt
	grep ">" $i |sed -e 's/^>//' |cat - $Data/5_remove_Chimera_Greengene/$ID1  |cut -f1 |cut -d' ' -f1 > $Data/5_remove_Chimera_Greengene/$ID2
done

