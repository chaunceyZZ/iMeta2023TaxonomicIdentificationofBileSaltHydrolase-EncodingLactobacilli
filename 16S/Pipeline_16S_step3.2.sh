Data=$(dirname $PWD)
cd $Data;

for i in $(ls $Data/3_filtered_fa/*.fna);
do
	cat $i | grep -E ">" | wc -l;
done
