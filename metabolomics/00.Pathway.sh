#Analysis StepI: Extract Significant data 

#Data=/Users/caiyuanyuan/Desktop/FY/;
Data=/home/Project/Metabolome/LJK/LJK_Dis
cd $Data;

### Extract signficant data between groups using U test or KW test
[ ! -d $Data/4.10_ID_KEGG ] && mkdir $Data/4.10_ID_KEGG;
for i in $(ls $Data/ID_Metabolites/*.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
#        Rscript $Data/bin/A4.15_pathway_Enrich_analysis.r $i  $Data/bin/$ff-grouping.info $Data/4.10_ID_KEGG
	
done



[ ! -d $Data/4.11_ID_KEGG_Protein ] && mkdir $Data/4.11_ID_KEGG_Protein;
for i in $(ls $Data/ID_Metabolites/*.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        Rscript $Data/bin/A4.14_KEGGEnrich_analysis.r $i $Data/4.1_VIP_Sig/$ff.Utest.Sig.VIP.txt $Data/bin/$ff-grouping.info $Data/4.11_ID_KEGG_Protein
	
done
