use strict;
use File::Basename;
die("Argument: OTU_count_txt(otu.count.txt)\n") if ((@ARGV < 1) or (@ARGV >2));

my $rawFile=shift;
my ($fn, $dir, undef) = fileparse($rawFile, qw/.txt/);
$fn=~/(\w+)\.otu\.counts/;
my $out=$dir."/".$1."."."otu.abundance.txt";

open In,"<",$rawFile;
open Out,">",$out;


my %Data;
my %Indi;
my %Sum;
my @id;
while(<In>)
{
	chomp;
	if(/^ID/)
	{
		@id=split /\t/;
		for my $i(1..$#id)
		{
			$Indi{$id[$i]}=1;
		}
	}
	else
	{
		my @temp=split /\t/;
		for my $i(1..$#temp)
		{
			$Data{$temp[0]}{$id[$i]}=$temp[$i];
			$Sum{$id[$i]}=$Sum{$id[$i]}+$temp[$i];
		}
	}
}

printf Out "ID\t";
for my $i(sort keys %Indi)
{
	printf Out "%s\t",$i;
}

printf Out "\n";
for my $i(sort keys %Data)
{
	printf Out "%s\t",$i;
	for my $j(sort keys %Indi)
	{
		printf Out "%1.6f\t",$Data{$i}{$j}/$Sum{$j};
	}
	printf Out "\n";
}


