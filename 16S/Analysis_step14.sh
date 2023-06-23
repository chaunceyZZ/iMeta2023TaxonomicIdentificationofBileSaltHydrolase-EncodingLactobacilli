#Current working taxanomic level is Genus

#Data=/home/wmj/Metagenome-2/16s/ZSL-1
Data=$(dirname $PWD)

[ ! -d $Data/Analysis/18_raw_beta_PCoA ] && mkdir $Data/Analysis/18_raw_beta_PCoA;



for i in $(ls $Data/14_Statistics/beta/*/*.otu_table*.txt);
do
	dir1=${i#*beta}
	dir2=`echo "$dir1" | cut -d / -f 2`
	[ ! -d $Data/Analysis/18_raw_beta_PCoA/$dir2 ] && mkdir -p $Data/Analysis/18_raw_beta_PCoA/$dir2;
	cp $i $Data/Analysis/18_raw_beta_PCoA/$dir2;
	Rscript $Data/bin/A14.1_PCoA.r  $i  $Data/bin/$dir2-color.txt $Data/Analysis/18_raw_beta_PCoA/$dir2 2;		
done


for i in $(ls $Data/13_betaDiv/*/*.otu_table*.txt);
do
	b='pc.txt'
	if [[ $i == *$b$ ]]
	then
		echo "pc"
	else
		dir0=$(dirname $i)
		dir1=${dir0#*/13_betaDiv/}
		cp $i $Data/Analysis/18_raw_beta_PCoA/$dir1;
		Rscript $Data/bin/A14.2_PCoA.r  $i  $Data/bin/$dir2-color.txt $Data/Analysis/18_raw_beta_PCoA/$dir2 2;	
	fi	

done


