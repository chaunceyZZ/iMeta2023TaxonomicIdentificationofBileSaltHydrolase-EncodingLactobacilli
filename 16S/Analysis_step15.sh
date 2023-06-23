#Current working taxanomic level is Genus

Data=$(dirname $PWD)

[ ! -d $Data/Analysis/19_Sig_beta_PCoA3 ] && mkdir $Data/Analysis/19_Sig_beta_PCoA3;


for i in $(ls $Data/Analysis/12_Sig_Distance2/*/*.txt)
do
	dir1=$(basename $i)
	dir2=`echo "$dir1" | cut -d . -f 1`
	[ ! -d $Data/Analysis/19_Sig_beta_PCoA3/$dir2 ] && mkdir $Data/Analysis/19_Sig_beta_PCoA3/$dir2;
	Rscript $Data/bin/A15.1_Sig_PCoA.r  $i  $Data/Analysis/$dir2-color3.txt $Data/Analysis/19_Sig_beta_PCoA3/$dir2;	
done

