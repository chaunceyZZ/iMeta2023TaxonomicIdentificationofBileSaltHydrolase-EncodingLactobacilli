use strict;

use File::Basename;
die("Argument: InFile_1 InFile_2 OutDir\n") if ((@ARGV < 3) or (@ARGV >3));

my $InFile1=$ARGV[0];
my $InFile2=$ARGV[1];
my $dir=$ARGV[2];


my ($fn, $dir1, undef) = fileparse($InFile1, qw/.enriched.metabolite.dat/);

#my $dir="./";

open Ref1,"<",$InFile1;
my %Metabolite;
while(<Ref1>)
{
	chomp;
	if(/^M\d+/)
	{
		my @temp=split /\s+/;
		$Metabolite{$temp[0]}=$temp[1];
	}
}

open In1,"<",$InFile2;
open Out1,">",$dir."/".$fn.".enriched.gene.dat";
while(<In1>)
{
	my @temp=split /\s+/;
	printf Out1 "%s\t%s",$Metabolite{$temp[0]},$_;
}

close In1;
close Out1;

open In2,"<",$dir."/".$fn.".enriched.gene.dat";
open Out2,">",,$dir."/".$fn.".enriched.gene.nred.1.dat";
open Out3,">",,$dir."/".$fn.".enriched.gene.nred.2.dat";

my %NMeta;
my %NGene;
while(<In2>)
{
	chomp;
	my @temp=split /\s+/;
	if(exists $NMeta{$temp[1]}{$temp[5]})
	{
		if (!($NMeta{$temp[1]}{$temp[5]} eq $temp[0]))
		{
			$NMeta{$temp[1]}{$temp[5]}="NA";
		}
	}
	else
	{
		$NMeta{$temp[1]}{$temp[5]} = $temp[0];
	}

	if(exists $NGene{$temp[5]})
	{
		if(!($NGene{$temp[5]}))
		{
			$NGene{$temp[5]}="NA";
		}
	}
	else
	{
		$NGene{$temp[5]}=$temp[0];
	}
}

for my $m(sort keys %NMeta)
{
	my $count=1;
	for my $n(sort keys %{$NMeta{$m}})
	{
		if(!($NMeta{$m}{$n} eq "NA") && length($n)>0)
		{
			printf Out2 "%s\t%s\t%d\t%s\n",$NMeta{$m}{$n},$m,$count,$n;
			$count++;
		}
	}
}

for my $i(sort keys %NGene)
{
	if(!($NGene{$i} eq "NA") &&length($i)>0)
	{
		printf Out3 "%s\t%s\n",$NGene{$i},$i;
	}
}

close In2;
close Out2;
close Out3;
