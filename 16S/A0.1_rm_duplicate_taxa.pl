use strict;
use File::Basename;
die("Argument: tmp1_abundance.txt\n") if ((@ARGV < 1) or (@ARGV >2));

my $TaxaLevel="Genus";

my $rawFile=shift;
my ($fn, $dir, undef) = fileparse($rawFile, qw/.txt/);
$fn=~/(\w+)\.\w+\.abundance\.tmp1/;
my $outFile=$dir."/".$1.".".$TaxaLevel.".abundance.tmp2.txt";
open In,"<",$rawFile;
open Out,">",$outFile;

my %Data;
while(<In>)
{
	chomp;
	if(/^ID/)
	{
		printf Out "%s\n",$_;
	}
	else
	{
		/^(\S+)\s+(.*)/;
		my @temp=split /\s+/;
		my $total;
		for my $i(1..$#temp)
		{
			$total=$total+$temp[$i];
		}
		if(exists $Data{$1} && $Data{$1}{"abundance"}<$total)
		{
			$Data{$1}{"data"}=$2;
			$Data{$1}{"abundance"}=$total;
		}
		else
		{
			$Data{$1}{"data"}=$2;
			$Data{$1}{"abundance"}=$total;
		}
	}
}

for my $id(sort keys %Data)
{
	printf Out "%s\t%s\n",$id,$Data{$id}{"data"};
}




