#nalysis StepV: analysis for VIP and significant sites

Data=/home/Metagenome-2/Metabolome/HFQ/Pos/HFQ-Pos.P;
cd $Data;

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
	Rscript $Data/bin/A4.5_Plot_mds.r $i $Data/bin/$ff-grouping.info $Data/4.5_VIP_Sig_MDS 2
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
        perl $Data/bin/A4.8.1_MZ_HMDB_Gene.pl $i $Data/4.7_VIP_Sig_Anno;
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
