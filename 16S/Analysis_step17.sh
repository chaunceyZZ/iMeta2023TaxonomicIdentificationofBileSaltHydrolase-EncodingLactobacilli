#Current working taxanomic level is Genus

Data=$(dirname $PWD)

[ ! -d $Data/Analysis/21_raw_beta_PCA ] && mkdir $Data/Analysis/21_raw_beta_PCA;


# for i in $(ls $Data/14_Statistics/beta/*/*.otu_table*.txt);
# do
# 	dir1=${i#*beta}
# 	dir2=`echo "$dir1" | cut -d / -f 2`
# 	[ ! -d $Data/Analysis/21_raw_beta_PCA/$dir2 ] && mkdir -p $Data/Analysis/21_raw_beta_PCA/$dir2;
# 	cp $i $Data/Analysis/21_raw_beta_PCA/$dir2;
# 	Rscript $Data/bin/A18.1_PCA.r  $i  $Data/bin/$dir2-color.txt $Data/Analysis/21_raw_beta_PCA/$dir2 2;		
# done
# 
# 
# for i in $(ls $Data/13_betaDiv/*/*.otu_table*.txt);
# do
# 	b='pc.txt'
# 	if [[ $i == *$b$ ]]
# 	then
# 		echo "pc"
# 	else
# 		dir0=$(dirname $i)
# 		dir1=${dir0#*/13_betaDiv/}
# 		cp $i $Data/Analysis/21_raw_beta_PCA/$dir1;
# 		Rscript $Data/bin/A18.2_PCA.r  $i  $Data/bin/$dir2-color.txt $Data/Analysis/21_raw_beta_PCA/$dir2 2;	
# 	fi	
# 
# done
# 
[ ! -d $Data/Analysis/22_Sig_beta_PCA2 ] && mkdir $Data/Analysis/22_Sig_beta_PCA2;

for i in $(ls $Data/Analysis/12_Sig_Distance2/*/*.txt)
do
	dir1=$(basename $i)
	dir2=`echo "$dir1" | cut -d . -f 1`
	[ ! -d $Data/Analysis/22_Sig_beta_PCA2/$dir2 ] && mkdir $Data/Analysis/22_Sig_beta_PCA2/$dir2;
	Rscript $Data/bin/A18.1_Sig_PCA.r  $i  $Data/Analysis/$dir2-color2.txt $Data/Analysis/22_Sig_beta_PCA2/$dir2;	
done

