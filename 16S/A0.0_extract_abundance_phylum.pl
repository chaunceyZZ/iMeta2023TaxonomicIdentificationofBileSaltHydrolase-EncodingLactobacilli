use strict;
use File::Basename;
die("Argument: Original_genus_abundance(otu_table_L6.txt)\n") if ((@ARGV < 1) or (@ARGV >2));

my $TaxaLevel="Phylum";

my $rawFile=shift;
my ($fn, $dir, undef) = fileparse($rawFile, qw/.txt/);
$fn=~/(\w+)\.gut/;
my $outFile=$dir."/".$1.".".$TaxaLevel.".abundance.txt";
open In,"<",$rawFile;
open Out,">",$outFile;

while(<In>)
{
	chomp;
	if(/^\#OTU\s+(ID.*)/)
	{
		printf Out "%s\n",$1;
	}
	elsif(/^k_/)
	{
		my $line=$_;
		my @taxa=split /;/,(split /\s+/,$_)[0];	
		if(!($taxa[1] eq "Other") && !($taxa[1]=~/_$/) )
		{
#			$line=~/^\S+_(\S+)\s+(.*)/;
#			printf Out "%s\t%s\n",$1,$2;
#g__Candidatus Arthromitus 0
			$line=~/^\S+_(\S*\s??[A-Za-z]*)\s+(\d.*)/;
			my @phy=split /__/,(split /\s+/,$taxa[1])[0];
			printf Out "%s\t%s\n",$phy[1],$2;
			
		}

		
		
			
		
	}
}
