#Current working taxanomic level is Genus

Data=$(dirname $PWD)

##	MA/MAS/MC/MeH/MeS
[ ! -d $Data/Analysis/15_RDA_and_Heatmap ] && mkdir -p $Data/Analysis/15_RDA_and_Heatmap

for i in $(ls $Data/Analysis/4_common_OTUs/Total/*_table.txt)
do
	data=$(basename $i);
	
	Out_dir=$(dirname $i)
	d1=${Out_dir/4_common_OTUs/15_RDA_and_Heatmap};
	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	[ ! -d $d1 ] && mkdir -p $d1

	cp $i $d1
	sed -i "s/#/""/g" $d1/$data ### 删除"#"
	sed -i '1d' $d1/$data	###删除第一行



	arg1=("L.animalis-Control" "L.fermentum-Control" "L.salvarius-Control" "L.reuteri-Control" "L.reuteriAPS-Control" "APS-Control" )	####分组信息
	Rscript $Data/bin/A11.1_RDA.r $data $d1 ${arg1[*]}  $Data/bin/Total-color.txt;
	echo "A11.1"
	##参数为读入文件;保存路径;分组信息
	###### 主要生成p值
	p=0.001

	arg2=("L.animalis-Control" "L.fermentum-Control" "L.salvarius-Control" "L.reuteri-Control" "L.reuteriAPS-Control" "APS-Control" )  #####
	Rscript $Data/bin/A11.2_Info_P.r $d1 ${arg2[*]} $p  ${#arg1[*]} $Data/bin/Total-color.txt;	
	#参数为保存路径;分组信息;p值;分组数;原始的分组信息
	####### 挑选p值,全部小于P 
	echo "A11.2"

	arg3=("L.animalis-Control" "L.fermentum-Control" "L.salvarius-Control" "L.reuteri-Control" "L.reuteriAPS-Control" "APS-Control" )  ##### A-B:B高为UP-1,B低为down-0,一般B为model
	Rscript $Data/bin/A11.3_Info_Factor.r $d1 ${arg2[*]} ${#arg1[*]} $Data/bin/Total-color.txt;	
	#	factor:严格要求,同增同减
	# up 和down 以及 生成物种信息
	echo "A11.3"


	arg4=("Control" "L.animalis" "L.fermentum" "L.salvarius" "L.reuteri" "L.reuteriAPS" "APS")  ####Heatmap 
	###### 不画grid
	Rscript $Data/bin/A11.4_Heatmap_no_grid.r $d1 ${arg4[*]}  $Data/bin/Total-color.txt;
	### 保存路径;分组信息;具体组
	##### 画 grid

	echo "A11.4"
	arg5=("T1Control"	"T1Classical"	"T1T1H") #Grid;第一位必须是基准即疾病;第二位必须是正常组，后面是给药组
#	Rscript $Data/bin/A11.4_Heatmap_grid.r $d1 ${arg4[*]}  $Data/bin/Total-color.txt ${arg5[*]} $p;
	### 保存路径;分组信息;具体组;grid信息




done



