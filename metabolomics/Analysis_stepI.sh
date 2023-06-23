#Analysis StepI: Extract Significant data form LEfSe analysis
Data=/home/Metagenome-2/Metabolome/LJK/LJK_Pos_Dis;
#Data=/Users/caiyuanyuan/Desktop/FY/;
cd $Data;

###Plot Volcano figures
[ ! -d $Data/1.1_raw_Volcano ] && mkdir $Data/1.1_raw_Volcano;
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
	Rscript $Data/bin/A1.3_Plot_mds.r $i $Data/bin/$ff-grouping.info $Data/1.4_raw_MDS 2
	
done
echo "1.4_raw_MDS done"

###Get Heamap figure
[ ! -d $Data/1.5_raw_Heatmap ] && mkdir $Data/1.5_raw_Heatmap;
for i in $(ls $Data/0.0_raw_data/*.txt)
do
	ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
	Rscript $Data/bin/A1.4_Plot_heatmap.r $i $Data/bin/$ff-grouping.info $Data/1.5_raw_Heatmap 2;
	
done
echo "1.5_raw_Heatmap done"
###Convert mz data to HMDB_ID, to Gene_Name
#[ ! -d $Data/1.6_All_Anno ] && mkdir $Data/1.6_All_Anno;
#for i in $(ls $Data/0.0_raw_data/*.mz.dat)
#do
#        ID=$(basename $i);
#        perl $Data/bin/A1.5_MZ_HMDB_Gene.pl $i $Data/1.6_All_Anno;
#done

####Enrich analysis
#[ ! -d $Data/1.7_All_Enrich ] && mkdir $Data/1.7_All_Enrich;
#for i in $(ls $Data/0.0_raw_data/*.txt)
#do
#        ID=$(basename $i);
#        Rscript $Data/bin/A1.6_Group_enrich.r $i $Data/bin/${ID%.*.*}-grouping.info $Data/1.7_All_Enrich 2
#done

#for i in $(ls $Data/1.7_All_Enrich/*.enriched.metabolite.dat)
#do
#	ID=$(basename $i);
#	perl $Data/bin/A1.7_Enriched_metabolite_gene.pl $i $Data/1.6_All_Anno/${ID%.*.*.*}.gene.dat $Data/1.7_All_Enrich/;
#done

#for i in $(ls $Data/1.7_All_Enrich/*.enriched.gene.nred.2.dat)
#do
#        ID=$(basename $i);
#        Rscript $Data/bin/A1.8_Enrich_analysis.r $i $Data/1.7_All_Enrich/;
#done
