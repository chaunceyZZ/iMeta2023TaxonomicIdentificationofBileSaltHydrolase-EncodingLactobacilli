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


Depth_alpha=40000
Depth_beta=40000
############

########### NOTE
gr=("Total")
###########



for m in unifrac unweighted_unifrac weighted_unifrac;
do 
	for i in ${gr[@]}
	do
		make_2d_plots.py -i $Data/13_betaDiv/$i/${m}_$i.gut.sub3W1.filter_4_noChimera_rmAlignFail.otu_table_even${Depth_beta}.pc.txt -m $Data/bin/$i-grouping.info -o $Data/13_betaDiv/$i/ ;
	done
done

