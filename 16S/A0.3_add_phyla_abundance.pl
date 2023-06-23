use strict;
use File::Basename;
die("Argument: Original_genus_abundance(otu_table_L6.txt)\n") if ((@ARGV < 1) or (@ARGV >2));

my $refFile=shift;
my ($fn, $dir, undef) = fileparse($refFile, qw/.txt/);

open Ref,"<",$refFile;
my %data;
while(<Ref>)
{
        if(!/^#/)
        {
                my @temp=split /\s+/;
                if($temp[0]=~/p__(\S+)\;c__.*g__(\S+)/)
                {
                        $data{$2}=$1;
                }
        }
}

$fn=~/(\w+)\.gut/;
open In,"<",$dir."/".$1.".Genus.abundance.txt";
open Out,">",$dir."/".$1.".PlyGen.abundance.txt";
my %NewData;
while(<In>)
{
	my $line=$_;
	if(/^ID/)
	{
		$line=~s/ID/Genera/;
		printf Out "Phyla %s",$line;
	}
	else
	{
		my @temp=split /\s+/;
		$NewData{$data{$temp[0]}}{$temp[0]}=$line;
	}
}
for my $ply(sort keys %NewData)
{
	for my $gen(sort keys %{$NewData{$ply}})
	{
		printf Out "%s %s",$ply,$NewData{$ply}{$gen};
	}
}
