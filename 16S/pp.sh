#StepIV:Alignment and remove Chimera
#NOTE: The QIIME enviroment should be set firstly!


Data=/home/Project/16S/20200225_WXD
GREENGENE=/home/Project/Ref_database/greengenes;


Depth_alpha=23000
Depth_beta=23000
############

########### NOTE
gr=("Total")
###########

[ ! -d $Data/12_alphaRarefactionPlot ] && mkdir $Data/12_alphaRarefactionPlot
for i in ${gr[@]}
do
	make_rarefaction_plots.py -i $Data/11_collated_alpha/$i -m $Data/bin/$i-grouping.info -g pdf -d 180 -o $Data/12_alphaRarefactionPlot/$i --generate_average_tables
done
