Data=/home/Metagenome-2/Metabolome/LJK/LJK_Pos_Dis;
#Data=/Users/caiyuanyuan/Desktop/FY/;
cd $Data;


##	ID_Metabolites ¼ø¶¨Àë×Ó
# 5.0_ID_Metabolites

### Extract signficant data between groups using U test or KW test
[ ! -d $Data/5.1_ID_KEGG ] && mkdir $Data/5.1_ID_KEGG;
for i in $(ls $Data/5.0_ID_Metabolites/*.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        Rscript $Data/bin/A4.15_pathway_Enrich_analysis.r $i  $Data/bin/$ff-grouping.info $Data/5.1_ID_KEGG
	
done

[ ! -d $Data/5.2_ID_KEGG_Protein ] && mkdir $Data/5.2_ID_KEGG_Protein;
for i in $(ls $Data/5.0_ID_Metabolites/*.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        Rscript $Data/bin/A4.14_KEGGEnrich_analysis.r $i $Data/4.1_VIP_Sig/$ff.Utest.Sig.VIP.txt $Data/bin/$ff-grouping.info $Data/5.2_ID_KEGG_Protein
	
done



