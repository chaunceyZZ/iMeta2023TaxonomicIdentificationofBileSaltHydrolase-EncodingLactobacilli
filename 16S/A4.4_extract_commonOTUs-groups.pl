use strict;
use File::Basename;
#die("Argument: OTU_Abundance|Count_File(otu.abundance|counts.txt) Grouping_File(*-grouping.info) Column\n") if ((@ARGV < 3) or (@ARGV >4));

my $InFile=$ARGV[0];
my $GroupFile=$ARGV[1];
my $grpCol=$ARGV[2]-1;

my ($fn, $dir, undef) = fileparse($InFile, qw/.txt/);
$fn=~/(\w+)\.otu\.(\w+)/;

my $L1=$1;
my $L2=$2;

my %Groups;
my $level;
open Ref,"<",$GroupFile;
while(<Ref>)
{
	if(/ID/)
	{
		my @temp=split /\s+/;
		$level=$temp[$grpCol];
	}
	else
	{
		my @temp=split /\s+/;
		$Groups{$temp[$grpCol]}{$temp[0]}=1;
	}
}

for (my $cutoff=0.05;$cutoff<=1;$cutoff=$cutoff+0.05)
{
	open In,"<",$InFile;
	open Out,">",$dir.$L1."/".$L1.".commonOTUs-".$level."-".$cutoff.".".$L2.".txt";

	my @Indi;

	while(<In>)
	{
		chomp;
		if(/^ID/)
		{
			@Indi=split /\s+/;
			printf Out "%s\n",$_;
		}
		else
		{
			my $line=$_;
			my @temp=split /\s+/;

			my %Data;
			for my $i(1..$#temp)
			{
				$Data{$Indi[$i]}=$temp[$i];
			}
			my $flag=1;
			for my $Gid(sort keys %Groups)
			{	
				my ($tmp_counts,$tmp_total)=(0,0);
				for my $Iid(sort keys %{$Groups{$Gid}})
				{
					if($Data{$Iid} > 0)
					{
						$tmp_counts++;
					}
					$tmp_total++;
				}
				if($tmp_counts/$tmp_total < $cutoff)
				{
					$flag=0;
				}
			}
			if($flag == 1 )
			{
				printf Out "%s\n",$line;
			}
		}
	}
}

