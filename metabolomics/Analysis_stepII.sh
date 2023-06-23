#Analysis StepI: Extract Significant data 

#Data=/Users/caiyuanyuan/Desktop/FY/;
Data=/home/Metagenome-2/Metabolome/HFQ/Pos/HFQ-Pos.P;

cd $Data;


### Extract signficant data between groups using U test or KW test
[ ! -d $Data/2.1_Sig_KWU ] && mkdir $Data/2.1_Sig_KWU;
for i in $(ls $Data/0.0_raw_data/*.data.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        Rscript $Data/bin/A2.1_Extract_Sig_KWU.r $i $Data/bin/$ff-grouping.info $Data/2.1_Sig_KWU 2;
	
done

for i in $(ls $Data/2.1_Sig_KWU/*Sig.data.txt);
do
        perl -pi -e "s/^\t/ID\t/g" $i;
done

###Plot Volcano figure
[ ! -d $Data/2.2_Sig_Volcano ] && mkdir $Data/2.2_Sig_Volcano;
for i in $(ls $Data/2.1_Sig_KWU/*.Utest.Sig.data.txt)
do
	ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
	if [ "$ff" = "Fall" ]
	then
		continue
	else
		Rscript $Data/bin/A2.2_Plot_Volcano.r $i $Data/bin/$ff-grouping.info $Data/2.2_Sig_Volcano 2 ;
	fi
	echo "2.2_Sig_Volcano done"
done
#done

### Prepare data for LEfSe
[ ! -d $Data/2.3_Sig_LDA ] && mkdir $Data/2.3_Sig_LDA;
for i in $(ls $Data/0.0_raw_data/*.data.txt)
do
	ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
	perl $Data/bin/A2.3_Prepare_LEfSe.pl $i $Data/2.3_Sig_LDA/$ff.4lefse.txt $Data/bin/$ff-grouping.info 2;
done

###	Run LEfSe
for i in $(ls $Data/2.3_Sig_LDA/*.4lefse.txt);
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        format_input.py $i $Data/2.3_Sig_LDA/$ff.in -c 1 -u 2 -o 1000000;
        run_lefse.py $Data/2.3_Sig_LDA/$ff.in $Data/2.3_Sig_LDA/$ff.res;
        plot_res.py $Data/2.3_Sig_LDA/$ff.res --format pdf $Data/2.3_Sig_LDA/$ff.pdf;
        plot_features.py  --format pdf --archive zip  $Data/2.3_Sig_LDA/$ff.in $Data/2.3_Sig_LDA/$ff.res $Data/2.3_Sig_LDA/$ff.zip;
done

#Extract signficant data by LDA score
#for i in $(ls $Data/2.3_Sig_LDA/*.res);
#do
#        ID=$(basename $i);
#        perl $Data/bin/A2.4_Extract_Sig_LDA.pl $Data/0.0_raw_data/${ID%.*}.data.txt $i $Data/2.3_Sig_LDA/${ID%.*}.LDA.Sig.data.txt;
#done

### Calculate Distance for Significant Data
[ ! -d $Data/2.4_Sig_Distance ] && mkdir $Data/2.4_Sig_Distance;
for i in $(ls $Data/2.1_Sig_KWU/*Sig.data.txt)
do
       ID=$(basename $i);
       Rscript $Data/bin/A2.5_Calculate_Distance.r $i $Data/2.4_Sig_Distance;
done

###Get PCA Plot
[ ! -d $Data/2.5_Sig_PCA ] && mkdir $Data/2.5_Sig_PCA;
for i in $(ls $Data/2.4_Sig_Distance/*.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        Rscript $Data/bin/A2.6_Plot_pca.r $i $Data/bin/$ff-grouping.info $Data/2.5_Sig_PCA 2
done

###Get MDS Plot
[ ! -d $Data/2.6_Sig_MDS ] && mkdir $Data/2.6_Sig_MDS;
for i in $(ls $Data/2.4_Sig_Distance/*.txt)
do
	ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
	Rscript $Data/bin/A2.7_Plot_mds.r $i $Data/bin/$ff-grouping.info $Data/2.6_Sig_MDS 2
done

###Get Heatmap figure

[ ! -d $Data/2.7_Sig_Heatmap ] && mkdir $Data/2.7_Sig_Heatmap;
for i in $(ls $Data/2*/*Sig.data.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        Rscript $Data/bin/A2.8_Plot_heatmap.r $i $Data/bin/$ff-grouping.info $Data/2.7_Sig_Heatmap 2
done

#Convert mz data to HMDB_ID, to Gene_Name
[ ! -d $Data/2.8_Sig_Anno ] && mkdir $Data/2.8_Sig_Anno;
for i in $(ls $Data/2*/*Sig.data.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        perl $Data/bin/A2.9_Extract_SIG_MZ.pl $i $Data/0.0_raw_data/$ff.mz.dat $Data/2.8_Sig_Anno/$ff.mz.dat;
done
for i in $(ls $Data/2.8_Sig_Anno/*.mz.dat)
do
        ID=$(basename $i);
        perl $Data/bin/A2.10_MZ_HMDB_Gene.pl $i $Data/2.8_Sig_Anno;
done


####Enrich analysis
[ ! -d $Data/2.9_Sig_Enrich ] && mkdir $Data/2.9_Sig_Enrich;
for i in $(ls $Data/2*/*Sig.data.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        Rscript $Data/bin/A2.11_Group_enrich.r $i $Data/bin/$ff-grouping.info $Data/2.9_Sig_Enrich 2
done

for i in $(ls $Data/2.9_Sig_Enrich/*.enriched.metabolite.dat)
do
	ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        perl $Data/bin/A2.12_Enriched_metabolite_gene.pl $i $Data/2.8_Sig_Anno/$ff.gene.dat $Data/2.9_Sig_Enrich/;
done

for i in $(ls $Data/2.9_Sig_Enrich/*.enriched.gene.nred.2.dat)
do
        ID=$(basename $i);
        Rscript $Data/bin/A2.13_Enrich_analysis.r $i $Data/2.9_Sig_Enrich/;
done

#####Metabolite Network
[ ! -d $Data/2.10_Sig_Metbo_NW ] && mkdir $Data/2.10_Sig_Metbo_NW;
for i in $(ls $Data/2.8_Sig_Anno/*.HMDB.dat)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        perl $Data/bin/A2.14_Extract_SIG_Connect.pl $i $Data/2.10_Sig_Metbo_NW/$ff.connect;
done

#####Metabolite Enriched Pathway
for i in $(ls $Data/2.10_Sig_Metbo_NW/*.pathsum)
do
	Rscript $Data/bin/A2.15_Metobolite_enrich.r $i $Data/2.10_Sig_Metbo_NW/;
done     ID=$(basename $i);


