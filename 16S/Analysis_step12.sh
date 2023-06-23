#Current working taxanomic level is Genus

Data=$(dirname $PWD)


#cp $Data/12_alphaRarefactionPlot/averSSage_tables/*Group.txt $Data/Analysis/16_alphaRarefactionPlot/0_pre

####	运行R
## 曲线
gr=("Total")
for L in ${gr[@]}
do
for i in $(ls $Data/12_alphaRarefactionPlot/$L/average_tables/*Group.txt);
do
	
	d2=$Data/Analysis/16_alphaRarefactionPlot/$L
	echo $d2
	[ ! -d $d2 ] && mkdir -p $d2;
	cp $i $d2
	file1=$(basename $i)

	sed -i '/#/d' $d2/$file1 ##删除含#的文件，并直接作用于源文件（-i）
	sed -i '/color/d' $d2/$file1
#	sed -i '/xmax/d' $i
	sed -i 's/>> //g' $d2/$file1 #全局（g）替换（s）把前者替换为后者
	sed -i 's/xaxis:/xaxis 	/g' $d2/$file1
	sed -i 's/xmax:/xmax	/g' $d2/$file1
	sed -i 's/series/series	/g' $d2/$file1
	sed -i 's/error/error	/g' $d2/$file1

	Rscript $Data/bin/A12.1_alphaplot.r $d2  $d2/$file1 $Data/bin/$L-color.txt
	
	
done
done

####	运行R


## boxplot
for i in $(ls $Data/Analysis/7_OTUs_alpha/*/*.final.otu.*.txt)
do
	Out_dir=$(dirname $i)
	d1=${Out_dir/7_OTUs_alpha/16_alphaRarefactionPlot};
	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	[ ! -d $d1 ] && mkdir -p $d1;
	cp $i $d1
	Rscript $Data/bin/A12.2_alpha_boxplot.r $i $Data/bin/$c1-color.txt $d1
	
done






