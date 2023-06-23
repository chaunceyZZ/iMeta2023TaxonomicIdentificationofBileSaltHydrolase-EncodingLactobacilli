#Current working taxanomic level is Genus

Data=$(dirname $PWD)

[ ! -d $Data/Analysis/20_PICRUSt_Function ] && mkdir $Data/Analysis/20_PICRUSt_Function;

cp $Data/Analysis/10_LEfSe/*/*kegg.L3.res $Data/Analysis/20_PICRUSt_Function;


for i in $(ls $Data/Analysis/9_PICRUSt/*/KEGG_Total/*.kegg.L2_L2.txt)
do
	data=$(basename $i);
	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	[ ! -d $d1 ] && mkdir -p $d1

	cp $i $Data/Analysis/20_PICRUSt_Function;
	sed -i "s/#/""/g" $Data/Analysis/20_PICRUSt_Function/$data ### 删除"#"
	sed -i '1d' $Data/Analysis/20_PICRUSt_Function/$data	###删除第一行

	ff=$Data/Analysis/20_PICRUSt_Function/$c1".kegg.L3.res"

	Rscript $Data/bin/A17.1_Picrust_Function.r  $Data/Analysis/20_PICRUSt_Function/$data  $ff $Data/bin/$c1-color.txt $Data/Analysis/20_PICRUSt_Function;	

done



