#Analysis StepI: Extract Significant data 

#Data=/Users/caiyuanyuan/Desktop/FY/;
Data=/project/Metabolome/HFQ/HFQ_New-1/P_Dis
cd $Data;

[ ! -d $Data/4.11_ID_KEGG_Protein ] && mkdir $Data/4.11_ID_KEGG_Protein;
for i in $(ls $Data/ID_Metabolites/*.txt)
do
        ID=$(basename $i);
	ff=`echo "$ID" | cut -d . -f 1`
        Rscript $Data/bin/A4.14_KEGGEnrich_analysis.r $i $Data/4.1_VIP_Sig/$ff.Utest.Sig.VIP.txt $Data/bin/$ff-grouping.info $Data/4.11_ID_KEGG_Protein
	
done

