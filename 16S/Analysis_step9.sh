#Analysis StepI: Classify the abundance data 
#Current working taxanomic level is Genus

Data=$(dirname $PWD)
cd $Data;

[ ! -d $Data/Analysis/10_LEfSe ] && mkdir $Data/Analysis/10_LEfSe;


for i in $(ls $Data/Analysis/4_common_OTUs/*/*.otu.counts.txt)
do
	dd=$(dirname $i)
	d10=${dd/4_common_OTUs/10_LEfSe};

	ff=$(basename $i)
	c1=`echo "$ff" | cut -d . -f 1`


	[ ! -d $d10 ] && mkdir -p $d10;
	

	cp  $Data/Analysis/1_Abundance/*.Genus.abundance.txt $d10
	d9=${dd/4_common_OTUs/9_PICRUSt};
	cp $d9/*.L3.abundance.txt $d10


done


for i in $(ls $Data/Analysis/9_PICRUSt/*/*.normolized.out.biom);
do 
	ID=$(basename $i);
	dd=$(dirname $i)
	d10=${dd/9_PICRUSt/10_LEfSe};
	biom convert -i $i -o $d10/${ID%.*}.tab --table-type "OTU table" --to-tsv --header-key="taxonomy" --output-metadata-id="Lineage";
done

for i in $(ls $Data/Analysis/10_LEfSe/*/*.tab);
do
	perl $Data/bin/A9.0_Convert_OTU_table.pl $i;
done



#Prepare the abundance data for LEfSe
for i in $(ls $Data/Analysis/10_LEfSe/*/*.L3.abundance.txt);
do
	ID=$(basename $i);
	dd=$(dirname $i)
	ff=$(basename $i)
	c1=`echo "$ff" | cut -d . -f 1`
	perl $Data/bin/A9.1_Prepare_LEfSe.pl $i $dd/${ID%.*.*}.4lefse.txt $Data/bin/$c1-color.txt 2;
done
for i in $(ls $Data/Analysis/10_LEfSe/*/*.out.abundance.txt);
do
        ID=$(basename $i);
	dd=$(dirname $i)
	ff=$(basename $i)
	c1=`echo "$ff" | cut -d . -f 1`

        perl $Data/bin/A9.1_Prepare_LEfSe.pl $i $dd/${ID%.*.*}.4lefse.txt $Data/bin/$c1-color.txt 2;
done

for i in $(ls $Data/Analysis/10_LEfSe/*/*.Genus.abundance.txt);
do
        ID=$(basename $i);
	dd=$(dirname $i)
	ff=$(basename $i)
	c1=`echo "$ff" | cut -d . -f 1`
        perl $Data/bin/A9.1_Prepare_LEfSe.pl $i $dd/${ID%.*.*}.4lefse.txt $Data/bin/$c1-color.txt 2;
done

#Run LEfSe
for i in $(ls $Data/Analysis/10_LEfSe/*/*.4lefse.txt);
do
	ID=$(basename $i);
	dd=$(dirname $i)
	lefse-format_input.py $i $dd/${ID%.*.*}.in -c 1 -u 2 -o 1000000;
	run_lefse.py $dd/${ID%.*.*}.in $dd/${ID%.*.*}.res;
	lefse-plot_res.py $dd/${ID%.*.*}.res --format pdf $dd/${ID%.*.*}.pdf;
	lefse-plot_features.py --format pdf --archive zip $dd/${ID%.*.*}.in $dd/${ID%.*.*}.res $dd/${ID%.*.*}.zip;
done

for i in $(ls $Data/Analysis/10_LEfSe/*/*.out.4lefse.txt);
do
	ID=$(basename $i);
	dd=$(dirname $i)
	lefse-plot_cladogram.py $dd/${ID%.*.*}.res --format pdf $dd/${ID%.*.*}.clad.pdf --clade_sep 0.1;
	lefse-plot_cladogram.py $dd/${ID%.*.*}.res --format png $dd/${ID%.*.*}.clad.png --clade_sep 0.1;
done


######boxplot
for i in $(ls $Data/Analysis/9_PICRUSt/*/*.kegg.L3.abundance.txt);
do
	dd=$(dirname $i)
	d10=${dd/9_PICRUSt/10_LEfSe};
	ff=$(basename $i)
	c1=`echo "$ff" | cut -d . -f 1`
	Rscript $Data/bin/A9.2_Boxplot.r $i $d10 $dd/$c1.catagrized.kegg.L2.tab
done
