#nalysis StepI: Extract Significant data form LEfSe analysis

Data=/home/Metagenome-2/Metabolome/HFQ/Pos/HFQ-Pos.P;
cd $Data;

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

