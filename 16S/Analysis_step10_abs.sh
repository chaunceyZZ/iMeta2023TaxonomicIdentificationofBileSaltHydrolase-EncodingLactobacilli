#Current working taxanomic level is Genus

Data=$(dirname $PWD)
[ ! -d $Data/Analysis/11_Sig_abs ] && mkdir $Data/Analysis/11_Sig_abs;


### Extract signficant data between groups using U test or KW test

for i in $(ls $Data/Analysis/0_absolute_Abundance/*/Total.Genus.absolute.abundance.txt)
do

	Out_dir=$(dirname $i)
	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	d1=${Out_dir/0_absolute_Abundance/11_Sig_abs};
	[ ! -d $d1 ] && mkdir -p $d1;
	Rscript $Data/bin/A10.1_Extract_Sig_KWU_abs.r $i $Data/bin/$c1-color.txt $d1 2;
	#Data_File Group_File Out_Dir Group_column
done



for i in $(ls $Data/Analysis/11_Sig_abs/*/*Sig.data.txt);
do
       perl -pi -e "s/^\t/ID\t/g" $i;
done

### Calculate Distance for Significant Data
[ ! -d $Data/Analysis/12_Sig_Distance_abs ] && mkdir $Data/Analysis/12_Sig_Distance_abs;



for i in $(ls $Data/Analysis/11_Sig_abs/*/*Sig.data.txt)
do
	ID=$(basename $i);
	dd=$(dirname $i)
	d1=${dd/11_Sig_abs/12_Sig_Distance_abs}
	[ ! -d $d1 ] && mkdir -p $d1;
	Rscript $Data/bin/A10.2_Calculate_Distance.r $i $d1;

done

# ###Get PCA Plot
[ ! -d $Data/Analysis/13_Sig_PCA_abs ] && mkdir $Data/Analysis/13_Sig_PCA_abs;

for i in $(ls $Data/Analysis/12_Sig_Distance_abs/*/*.txt)
do

	dd=$(dirname $i)
       	d1=${dd/12_Sig_Distance_abs/13_Sig_PCA_abs}
	[ ! -d $d1 ] && mkdir -p $d1;
	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	# d2=${dd/12_Sig_Distance/4_common_OTUs}
        Rscript $Data/bin/A10.3_Plot_pca.r $i $Data/bin/$c1-color.txt $d1 2
done

###Get MDS Plot
[ ! -d $Data/Analysis/14_Sig_MDS_abs ] && mkdir $Data/Analysis/14_Sig_MDS_abs;

for i in $(ls $Data/Analysis/12_Sig_Distance_abs/*/*.txt)
do
       	dd=$(dirname $i)
       	d1=${dd/12_Sig_Distance_abs/14_Sig_MDS_abs}
	[ ! -d $d1 ] && mkdir $d1;
	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	# d2=${dd/12_Sig_Distance/4_common_OTUs}

        Rscript $Data/bin/A10.4_Plot_mds.r $i $Data/bin/$c1-color.txt $d1 2
done
