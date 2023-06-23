use strict;
use File::Basename;
die("Argument: OTU_count_table(otu_table.txt)\n") if ((@ARGV < 1) or (@ARGV >2));

my $rawFile=shift;
my ($fn, $dir, undef) = fileparse($rawFile, qw/.txt/);
$fn=~/(\w+)\.gut/;
my $out1=$dir."/".$1."."."otu.counts.txt";
my $out2=$dir."/".$1."."."OTU-Taxa.txt";

open In,"<",$rawFile;
open Out1,">",$out1;
open Out2,">",$out2;

while(<In>)
{
	chomp;
	if(/^\#OTU/)
	{
		my @temp=split /\t/;
		printf Out1 "ID\t";
		for my $i(1..($#temp-1))
		{
			printf Out1 "%s\t",$temp[$i];
		}
		printf Out1 "\n";
	}
	elsif(/^denovo/)
	{
		my @temp=split /\t/;
		for my $i(0..($#temp-1))
                {
                        printf Out1 "%s\t",$temp[$i];
                }
                printf Out1 "\n"; 
		printf Out2 "%s\t%s\n",$temp[0],$temp[$#temp];
	}
}
