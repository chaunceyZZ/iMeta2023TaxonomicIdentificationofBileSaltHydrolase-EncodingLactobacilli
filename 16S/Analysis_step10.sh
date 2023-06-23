#Current working taxanomic level is Genus

Data=$(dirname $PWD)
[ ! -d $Data/Analysis/11_Sig ] && mkdir $Data/Analysis/11_Sig;


### Extract signficant data between groups using U test or KW test

for i in $(ls $Data/Analysis/4_common_OTUs/*/*.abundance.txt)
do

	Out_dir=$(dirname $i)
	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	d1=${Out_dir/4_common_OTUs/11_Sig};
	[ ! -d $d1 ] && mkdir -p $d1;
	Rscript $Data/bin/A10.1_Extract_Sig_KWU.r $i $Data/bin/$c1-color.txt $d1 2;
	#Data_File Group_File Out_Dir Group_column
done



for i in $(ls $Data/Analysis/11_Sig/*/*Sig.data.txt);
do
       perl -pi -e "s/^\t/ID\t/g" $i;
done

### Calculate Distance for Significant Data
[ ! -d $Data/Analysis/12_Sig_Distance ] && mkdir $Data/Analysis/12_Sig_Distance;



for i in $(ls $Data/Analysis/11_Sig/*/*Sig.data.txt)
do
	ID=$(basename $i);
	dd=$(dirname $i)
	d1=${dd/11_Sig/12_Sig_Distance}
	[ ! -d $d1 ] && mkdir -p $d1;
	Rscript $Data/bin/A10.2_Calculate_Distance.r $i $d1;

done

###Get PCA Plot
[ ! -d $Data/Analysis/13_Sig_PCA ] && mkdir $Data/Analysis/13_Sig_PCA;

for i in $(ls $Data/Analysis/12_Sig_Distance/*/*.txt)
do

	dd=$(dirname $i)
       	d1=${dd/12_Sig_Distance/13_Sig_PCA}
	[ ! -d $d1 ] && mkdir -p $d1;
	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	d2=${dd/12_Sig_Distance/4_common_OTUs}
        Rscript $Data/bin/A10.3_Plot_pca.r $i $Data/bin/$c1-color.txt $d1 2
done

###Get MDS Plot
[ ! -d $Data/Analysis/14_Sig_MDS ] && mkdir $Data/Analysis/14_Sig_MDS;

for i in $(ls $Data/Analysis/12_Sig_Distance/*/*.txt)
do
       	dd=$(dirname $i)
       	d1=${dd/12_Sig_Distance/14_Sig_MDS}
	[ ! -d $d1 ] && mkdir $d1;
	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	d2=${dd/12_Sig_Distance/4_common_OTUs}

       Rscript $Data/bin/A10.4_Plot_mds.r $i $Data/bin/$c1-color.txt $d1 2
done
