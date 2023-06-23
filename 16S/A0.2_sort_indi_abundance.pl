use strict;
use File::Basename;
die("Argument: tmp2_abundance.txt\n\n") if ((@ARGV < 1) or (@ARGV >2));

my $TaxaLevel="Genus";

my $rawFile=shift;
my ($fn, $dir, undef) = fileparse($rawFile, qw/.txt/);
$fn=~/(\w+)\.\w+\.abundance\.tmp2/;
my $outFile=$dir."/".$1.".".$TaxaLevel.".abundance.txt";
open In,"<",$rawFile;
open Out,">",$outFile;

my %Data;
my %Indi;
my @id;
while(<In>)
{
	chomp;
	if(/^ID/)
	{
		@id=split /\s+/;
		for my $i(1..$#id)
		{
			$Indi{$id[$i]}=1;
		}
	}
	else
	{
		my @temp=split /\s+/;
		for my $i(1..$#temp)
		{
			$Data{$temp[0]}{$id[$i]}=$temp[$i];
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
		printf Out "%1.5f\t",$Data{$i}{$j};
	}
	printf Out "\n";
}
