#StepIII:Pick OTUs
#NOTE: This step could be taken days!
#NOTE: The QIIME enviroment should be set firstly!
### 选择OTU
Data=$(dirname $PWD)
[ ! -d $Data/4_Picked_OTUs ] && mkdir $Data/4_Picked_OTUs;

####
###	sub3W1 :note modify
####


#3.1 For Gut data
for i in $(ls $Data/3_filtered_fa/*gut*)
do
	echo "1"
	pick_otus.py -i $i -o $Data/4_Picked_OTUs/;  ##从fasta文件中提取OTUs
	ID1=$(basename ${i%.fa})_otus.txt
	ID2=$(basename ${i%.fa})_rep_set.fa
	echo "2"
	pick_rep_set.py -i $Data/4_Picked_OTUs/$ID1 -f $i -o $Data/4_Picked_OTUs/$ID2; 
	###选取代表序列,由于每个OTU中的序列不完全相同，因此需要选取一条代表性序列作为该OTU的序列，用于后续分析
done
 
echo "3"
#3.2 Keep the OTUs with minimum 4 reads.
for i in $(ls $Data/4_Picked_OTUs/*.gut.sub3W1_otus.txt)
do
	ID2=$(basename ${i%_otus.txt})_rep_set.fa
	$Data/bin/filter_OTUs_byReads.pl $i $Data/4_Picked_OTUs/$ID2 4;
done


