#StepIV:Alignment and remove Chimera
#NOTE: The QIIME enviroment should be set firstly!
#### docker run -it -v `pwd`:/home --name=xxxx yoshikiv/basespace-qiime-191-dev
Data=$(dirname $PWD)
# GREENGENE=/home/Project/Ref_database/greengenes;
GREENGENE=$Data/bin/greengenes


####
###	sub3W1 :note modify
####


###
###	gut-grouping:#SampleID Group
###


###
###	gr=(Total,MC,MeS,MeH,MA,MAS) :note modify
###

############	Note:Less than SubsampleSeq,if too large,will be missing files in 12_alphaRarefactionPlot;less 1W
##### Note
###	Get min
#	for i in $(ls $Data/3_filtered_fa/*.fa);do cat $i | grep -E ">" | wc -l;done

###	for i in $(ls $Data/3_filtered_fa/*.fa);do cat $i | grep -E ">" | wc -l;done


Depth_alpha=23000
Depth_beta=23000
############

########### NOTE
gr=("Total")
###########



echo "12_alphaRarefactionPlot"
[ ! -d $Data/12_alphaRarefactionPlot ] && mkdir $Data/12_alphaRarefactionPlot
for i in ${gr[@]}
do
	make_rarefaction_plots.py -i $Data/11_collated_alpha/$i -m $Data/bin/$i-grouping.info -g pdf -d 180 -o $Data/12_alphaRarefactionPlot/$i --generate_average_tables
done
