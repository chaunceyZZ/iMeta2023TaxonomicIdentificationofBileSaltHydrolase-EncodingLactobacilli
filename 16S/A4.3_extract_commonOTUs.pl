use strict;
use File::Basename;
die("Argument: OTU_Abundance or Count File(otu.abundance|counts.txt)\n") if ((@ARGV < 1) or (@ARGV >2));

my $rawFile=shift;
my ($fn, $dir, undef) = fileparse($rawFile, qw/.txt/);
$fn=~/(\w+)\.otu\.(\w+)/;

my $L1=$1;
my $L2=$2;

for (my $cutoff=0.05;$cutoff<=1;$cutoff=$cutoff+0.05)
{
	open In,"<",$rawFile;
	open Out,">",$dir.$L1."/".$L1.".commonOTUs-".$cutoff.".".$L2.".txt";
	
	while(<In>)
	{
		chomp;
		if(/^ID/)
		{
			printf Out "%s\n",$_;
		}
		else
		{
			my $line=$_;
			my @temp=split /\s+/;
			my ($tmp_counts,$tmp_total)=(0,0);
			for my $i(1..$#temp)
			{
				if($temp[$i] >0)
				{
					$tmp_counts++;
				}
				$tmp_total++;
			}
			if($tmp_counts/$tmp_total > $cutoff)
			{
				printf Out "%s\n",$line;
			}
		}
	}
}

