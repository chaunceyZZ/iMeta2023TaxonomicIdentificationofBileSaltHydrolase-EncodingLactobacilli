use strict;
use File::Basename;
die("Argument: InFile RefFile OutFile\n") if ((@ARGV < 2) or (@ARGV >2));

my $InFile=$ARGV[0];
my $OutFile3=$ARGV[1];

my ($fn, $dir1, undef) = fileparse($OutFile3, qw/.connect/);


#my $InFile="cont.LDA.Sig.VIP.HMDB.dat";
#my $OutFile="cont.LDA.Sig.VIP.HMDB.connect";

my $OutFile1=$dir1."/".$fn.".path";
my $OutFile2=$dir1."/".$fn.".pathsum";

open Ref,"<","/home/Metagenome-2/Ref_database/HMDB/Metabolite_Pathway.txt";
my %Meta;
while(<Ref>)
{
	chomp;
	my @temp=split /\t/;
	$Meta{$temp[0]}{$temp[1]}=$temp[2];
}

open In,"<",$InFile;
open Out1,">",$OutFile1;
my %Data;
my %PathSum;
my %IndiMeta;
my $all;
while(<In>)
{
	chomp;
	my @temp=split /\s+/;
	if($temp[0]=~/^M\d+/ && $temp[2]=~/^HMDB/)
	{
		if(exists $Meta{$temp[2]})
		{
			for my $id(sort keys %{$Meta{$temp[2]}})
			{
				$Data{$temp[0]}{$temp[2]}{$id}=$Meta{$temp[2]}{$id};
				printf Out1 "%s\t%s\t%s\t%s\n",$temp[0],$temp[2],$id,$Meta{$temp[2]}{$id};

				$PathSum{$Meta{$temp[2]}{$id}}=$PathSum{$Meta{$temp[2]}{$id}} +1;
				$IndiMeta{$temp[0]}=1;
			}
		}
	}
}
$all=keys %IndiMeta;

open Out2,">",$OutFile2;
for my $id(sort keys %PathSum)
{
	printf Out2 "%s\t%d\t%d\n",$id,$PathSum{$id},$all;
}

my %Connect;
for my $Ion1(sort keys %Data)
{
	for my $Ion2(sort keys %Data)
	{
		if(!($Ion1 eq $Ion2))
		{
			for my $ID1(sort keys %{$Data{$Ion1}})
			{
				for my $P1(sort keys %{$Data{$Ion1}{$ID1}})
				{
					for my $ID2(sort keys %{$Data{$Ion2}})
					{
						for my $P2(sort keys %{$Data{$Ion2}{$ID2}})
						{
							if($Data{$Ion1}{$ID1}{$P1} eq $Data{$Ion2}{$ID2}{$P2})
							{
								$Connect{$Ion1}{$Ion2}{$ID1."_".$ID2}{$P1."_".$P2}=$Data{$Ion2}{$ID2}{$P2};							  
							}
						}
					}
				}
			}
		}
	}
}

open Out3,">",$OutFile3;
for my $Ion1(sort keys %Connect)
{
	for my $Ion2(sort keys %{$Connect{$Ion1}})
	{
		for my $id(sort keys %{$Connect{$Ion1}{$Ion2}})
		{
			for my $p(sort keys %{$Connect{$Ion1}{$Ion2}{$id}})
			{
				printf Out3 "%s\t%s\t%s\t%s\t%s\n",$Ion1,$Ion2,$id,$p,$Connect{$Ion1}{$Ion2}{$id}{$p};
			}
		}
	}
}

