Data=$(dirname $PWD)
cd $Data;
SIZE=400;
# SubsampleSeq=31000;
##############  NOTE:sub3W1

SubsampleSeq=`(for i in $(ls $Data/3_filtered_fa/*.fna);do cat $i | grep -E ">" | wc -l;done) | sort  -n | head -1`
echo $SubsampleSeq
