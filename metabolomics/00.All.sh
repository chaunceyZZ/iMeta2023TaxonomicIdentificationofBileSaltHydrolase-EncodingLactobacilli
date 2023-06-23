Data=/home/Project/TEMP/QFJ
#Data=/Users/caiyuanyuan/Desktop/FY/;
cd $Data;

###Plot Volcano figures
[ ! -d $Data/1.1_raw_Volcano ] && mkdir $Data/1.1_raw_Volcano;
#for i in $(ls $Data/0.0_raw_data/*.data.txt)
#do
#    ID=$(basename $i);
for i in $(ls $Data/0.0_raw_data/*.data.txt)
do
	ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
	if [ "$ff" = "Fall" ]
	then
		continue
	else
		Rscript $Data/bin/A1.0_Plot_Volcano.r $i $Data/bin/$ff-grouping.info $Data/1.1_raw_Volcano 2 ;
	fi
	
done
echo "1.1_raw_Volcano done"
#done

###Calculate various Distance
[ ! -d $Data/1.2_raw_Distance ] && mkdir $Data/1.2_raw_Distance;
for i in $(ls $Data/0.0_raw_data/*.data.txt)
do
	ID=$(basename $i);
	Rscript $Data/bin/A1.1_Calculate_Distance.r $i $Data/1.2_raw_Distance;
	
done
echo "1.2_raw_Distance done"

[ ! -d $Data/1.3_raw_PCA ] && mkdir $Data/1.3_raw_PCA;
for i in $(ls $Data/1.2_raw_Distance/*.txt)
do
	ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
    Rscript $Data/bin/A1.2_Plot_pca.r $i $Data/bin/$ff-grouping.info $Data/1.3_raw_PCA 2 $Data/bin/$ff-color.txt
	
done
echo "1.3_raw_PCA done"
###Get MDS Plot
[ ! -d $Data/1.4_raw_MDS ] && mkdir $Data/1.4_raw_MDS;
for i in $(ls $Data/1.2_raw_Distance/*.txt)
do
	ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
	Rscript $Data/bin/A1.3_Plot_mds.r $i $Data/bin/$ff-grouping.info $Data/1.4_raw_MDS 2 $Data/bin/$ff-color.txt
	
done
echo "1.4_raw_MDS done"

###Get Heamap figure
[ ! -d $Data/1.5_raw_Heatmap ] && mkdir $Data/1.5_raw_Heatmap;
for i in $(ls $Data/0.0_raw_data/*.txt)
do
	ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
	Rscript $Data/bin/A1.4_Plot_heatmap.r $i $Data/bin/$ff-grouping.info $Data/1.5_raw_Heatmap 2 $Data/bin/$ff-color.txt
	
done
echo "1.5_raw_Heatmap done"

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
#for i in $(ls $Data/2.1_Sig_KWU/*.data.txt)
#do
#        ID=$(basename $i);
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
        Rscript $Data/bin/A2.6_Plot_pca.r $i $Data/bin/$ff-grouping.info $Data/2.5_Sig_PCA 2 $Data/bin/$ff-color.txt
done

###Get MDS Plot
[ ! -d $Data/2.6_Sig_MDS ] && mkdir $Data/2.6_Sig_MDS;
for i in $(ls $Data/2.4_Sig_Distance/*.txt)
do
	ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
	Rscript $Data/bin/A2.7_Plot_mds.r $i $Data/bin/$ff-grouping.info $Data/2.6_Sig_MDS 2 $Data/bin/$ff-color.txt
done

###Get Heatmap figure

[ ! -d $Data/2.7_Sig_Heatmap ] && mkdir $Data/2.7_Sig_Heatmap;
for i in $(ls $Data/2*/*Sig.data.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        Rscript $Data/bin/A2.8_Plot_heatmap.r $i $Data/bin/$ff-grouping.info $Data/2.7_Sig_Heatmap 2 $Data/bin/$ff-color.txt
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
for i in $(ls $Data/2.1_Sig_KWU/*Sig.data.txt)
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
	ID=$(basename $i);
	Rscript $Data/bin/A2.15_Metobolite_enrich.r $i $Data/2.10_Sig_Metbo_NW/;
done   


### PLS-DA analysis
[ ! -d $Data/3.1_PLS_DA ] && mkdir $Data/3.1_PLS_DA;
for i in $(ls $Data/0.0_raw_data/*.data.txt)
do
	ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
	Rscript $Data/bin/A3.1_Plot_PLS_DA.r $i $Data/bin/$ff-grouping.info  $Data/3.1_PLS_DA 2 $Data/bin/$ff-color.txt;
done

### OPLS-DA by ropls
[ ! -d $Data/3.2_OPLS_DA ] && mkdir $Data/3.2_OPLS_DA;
cd $Data/3.2_OPLS_DA;
for i in $(ls $Data/0.0_raw_data/*.data.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        Rscript $Data/bin/A3.2_Plot_OPLS_DA.ropls.r $i $Data/bin/$ff-grouping.info  $Data/3.2_OPLS_DA 2 $Data/bin/$ff-color.txt;
        mv Rplots.pdf $Data/3.2_OPLS_DA/$ff.OPLSDA.summary.pdf;
        mv Rplots1.pdf $Data/3.2_OPLS_DA/$ff.OPLSDA.xscore.pdf;
	echo $ID "3.2_ropls"
done



### OPLS-DA by muma
cd $Data/3.2_OPLS_DA;
for i in $(ls $Data/0.0_raw_data/*.data.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        Rscript $Data/bin/A3.2_Plot_OPLS_DA.muma.r $i $Data/bin/$ff-grouping.info  $Data/3.2_OPLS_DA 2 $Data;
        mv $Data/3.2_OPLS_DA/OPLS-DAP/SPlot_OPLS-DA_P.pdf $Data/3.2_OPLS_DA/$ff.OPLSDA.splot.pdf;
        rm Rplot*.pdf;
        rm -r Groups/ OPLS-DAP/ PCA_Data_P/ Preprocessing_Data_P/;
	echo $ID "3.2_muma"
done



#Convert mz data to HMDB_ID, to Gene_Name

#Extract signficant data from VIP sites(OPLS-DA)
[ ! -d $Data/4.1_VIP_Sig ] && mkdir $Data/4.1_VIP_Sig ;
for i in $(ls $Data/2.1_Sig_KWU/*Sig.data.txt);
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        perl $Data/bin/A4.1_Extract_VIP_SIG.pl $i  $Data/3.2_OPLS_DA/$ff.OPLSDA.VIP.txt $Data/4.1_VIP_Sig/${ID%.*.*}.VIP.txt;
	echo $ID "4.1_Done"
done


###Plot Volcano figures
[ ! -d $Data/4.2_VIP_Sig_Volcano ] && mkdir $Data/4.2_VIP_Sig_Volcano;
for i in $(ls $Data/4.1_VIP_Sig/*.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
	if [ "$ff" = "Fall" ]
	then
		continue
	else
      		Rscript $Data/bin/A4.2_Plot_Volcano.r $i $Data/bin/$ff-grouping.info $Data/4.2_VIP_Sig_Volcano 2 ;
	fi
done

### Calculate Distance for VIP_Sig Data
[ ! -d $Data/4.3_VIP_Sig_Dis ] && mkdir $Data/4.3_VIP_Sig_Dis;
for i in $(ls $Data/4.1_VIP_Sig/*.txt)
do
       ID=$(basename $i);
       Rscript $Data/bin/A4.3_Calculate_Distance.r $i $Data/4.3_VIP_Sig_Dis;
	echo $ID "4.3_Done"
done

###Get PCA Plot
[ ! -d $Data/4.4_VIP_Sig_PCA ] && mkdir $Data/4.4_VIP_Sig_PCA;
for i in $(ls $Data/4.3_VIP_Sig_Dis/*.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        Rscript $Data/bin/A4.4_Plot_pca.r $i $Data/bin/$ff-grouping.info $Data/4.4_VIP_Sig_PCA 2 $Data/bin/$ff-color.txt

	echo $ID "4.4_Done"
done


###Get MDS Plot
[ ! -d $Data/4.5_VIP_Sig_MDS ] && mkdir $Data/4.5_VIP_Sig_MDS;
for i in $(ls $Data/4.3_VIP_Sig_Dis/*.txt)
do
	ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
	Rscript $Data/bin/A4.5_Plot_mds.r $i $Data/bin/$ff-grouping.info $Data/4.5_VIP_Sig_MDS 2 $Data/bin/$ff-color.txt
done

###Get Heatmap figure
[ ! -d $Data/4.6_VIP_Sig_Heatmap ] && mkdir $Data/4.6_VIP_Sig_Heatmap;
for i in $(ls $Data/4.1_VIP_Sig/*Sig.VIP.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`   
        Rscript $Data/bin/A4.6_Plot_heatmap.r $i $Data/bin/$ff-grouping.info $Data/4.6_VIP_Sig_Heatmap 2 $Data/bin/$ff-color.txt
done

#Convert mz data to HMDB_ID, to Gene_Name
[ ! -d $Data/4.7_VIP_Sig_Anno ] && mkdir $Data/4.7_VIP_Sig_Anno;
for i in $(ls $Data/4.1_VIP_Sig/*Sig.VIP.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        perl $Data/bin/A4.7_Extract_VIP_SIG_MZ.pl $i $Data/0.0_raw_data/$ff.mz.dat $Data/4.7_VIP_Sig_Anno/${ID%.*}.mz.dat;
	echo $ID "4.7_Done"

done

for i in $(ls $Data/4.7_VIP_Sig_Anno/*.mz.dat)
do
        ID=$(basename $i);
        perl $Data/bin/A4.8.1_MZ_HMDB_Gene.pl $i $Data/4.7_VIP_Sig_Anno $Data/bin;
done



####Enrich analysis
[ ! -d $Data/4.8_VIP_Sig_Enrich ] && mkdir $Data/4.8_VIP_Sig_Enrich;
cd $Data/4.8_VIP_Sig_Enrich;
for i in $(ls $Data/4.1_VIP_Sig/*Sig.VIP.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        Rscript $Data/bin/A4.9_Group_enrich.r $i $Data/bin/$ff-grouping.info $Data/4.8_VIP_Sig_Enrich 2
done

for i in $(ls $Data/4.8_VIP_Sig_Enrich/*.enriched.metabolite.dat)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        perl $Data/bin/A4.10_Enriched_metabolite_gene.pl $i $Data/4.7_VIP_Sig_Anno/$ff.VIP.gene.dat $Data/4.8_VIP_Sig_Enrich/;
done

for i in $(ls $Data/4.8_VIP_Sig_Enrich/*.enriched.gene.nred.2.dat)
do
 	ID=$(basename $i);
        Rscript $Data/bin/A4.11_Enrich_analysis.r $i $Data/4.8_VIP_Sig_Enrich/;
done
rm $Data/4.8_VIP_Sig_Enrich/*.xml

#####Metabolite Network
[ ! -d $Data/4.9_VIP_Sig_Metabo_NW ] && mkdir $Data/4.9_VIP_Sig_Metabo_NW;
for i in $(ls $Data/4.7_VIP_Sig_Anno/*.HMDB.dat)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        perl $Data/bin/A4.12_Extract_VIP_SIG_Connect.pl $i $Data/4.9_VIP_Sig_Metabo_NW/$ff.connect;
done

#####Metabolite Enriched Pathway
for i in $(ls $Data/4.9_VIP_Sig_Metabo_NW/*.pathsum)
do
        Rscript $Data/bin/A4.13_Metobolite_enrich.r $i $Data/4.9_VIP_Sig_Metabo_NW/;
done
