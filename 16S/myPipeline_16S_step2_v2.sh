#StepI:Convert the raw joined reads from Fastq format to Fasta format!
#NOTE: The QIIME enviroment should be set firstly!

### 将 Fastq格式转换成Fasta

Data=$(dirname $PWD)
cd $Data;

[ ! -d $Data/2_combine_fa ] && mkdir $Data/2_combine_fa;
cd $Data/2_combine_fa
:>commandlist_fastq2fasta
for i in $(ls $Data/1_combine_fq/*.fq);
do
	echo $(basename $i)
	echo "seqkit fq2fa --threads 4 $i -o $Data/2_combine_fa/$(basename $i .fq).fna;" >> commandlist_fastq2fasta
done
ParaFly -c commandlist_fastq2fasta -CPU 50


###convert_fastaqual_fastq.py
### 自带的质控
###可以将fastq拆分为fasta和quality; 也可以将将 fasta和quality合并为fastq
