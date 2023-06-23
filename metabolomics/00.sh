#nalysis StepI: Extract Significant data form LEfSe analysis
Data=/home/wmj/02_Metagenome/Metabolome/20190904_SZW/Temp
#Data=/Users/caiyuanyuan/Desktop/FY/;
cd $Data;

for i in $(ls $Data/4.7_VIP_Sig_Anno/*.mz.dat)
do
        ID=$(basename $i);
        perl $Data/bin/A4.8.1_MZ_HMDB_Gene.pl $i $Data/4.7_VIP_Sig_Anno;
done


