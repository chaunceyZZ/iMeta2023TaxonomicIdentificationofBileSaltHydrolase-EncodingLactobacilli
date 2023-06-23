use strict;

use File::Basename;
die("Argument: InFile OutDir\n") if ((@ARGV < 3) or (@ARGV >3));

#my $InFile="low.LDA.Sig.mz.dat";
#my $dir="./";

my $InFile=$ARGV[0];
my ($fn, $dir1, undef) = fileparse($InFile, qw/.mz.dat/);
my $dir=$ARGV[1];
my $dir0=$ARGV[2];


my %Data;
open Ref1,"<",$dir0."/megadb.20190908.txt";
while(<Ref1>)
{
		chomp;
        my @temp=split /\t/;
        $Data{$temp[0]}{"mz"}=$temp[3];
        $Data{$temp[0]}{"chem"}=$temp[4];
}
open Ref2,"<",$dir0."/HMDB.20190908.txt";
while(<Ref2>)
{		
		chomp;
        my @temp=split /\t/;
        $Data{$temp[0]}{"id"}=$temp[1];
		$Data{$temp[0]}{'Kingdom'}=$temp[2];
		$Data{$temp[0]}{'Super_Class'}=$temp[3];
		$Data{$temp[0]}{'Class'}=$temp[4];
		$Data{$temp[0]}{'Sub_Class'}=$temp[5];
		$Data{$temp[0]}{'Source'}=$temp[6];
		$Data{$temp[0]}{'pathway'}=$temp[11];
## Num	HMDB_ID	Kingdom	Super_Class	Class	Sub_Class	Source/Biological_role	Process	Cellular_Locations	Biospecimen_Locations	substituent	pathway	web

}
my %HMDB;
open Ref3,"<",$dir0."/ALL_HMDB_name.txt";
while(<Ref3>)
{
		chomp;
        my @temp=split /\t/;
        $HMDB{$temp[0]}=$temp[1];
}

open In1,"<",$InFile;
open Out1,">",$dir."/".$fn.".HMDB.dat";
printf Out1 "ID\tHit_num\tHMDB_ID\tChem\tName\tKingdom\tSuper_Class\tClass\tSub_Class\tSource\tPathway\tObs_mz\tRef_mz\n";

while(<In1>)
{
        chomp;
        if(! /name/)
        {
                my @temp=split /\t/;
                my $count=0;
                for my $id(sort keys %Data)
                {
                        if($Data{$id}{"mz"}>= $temp[2] && $Data{$id}{"mz"}<= $temp[3])
                        {
                                $count++;
                                printf Out1 "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$temp[0],$count,$Data{$id}{"id"},$Data{$id}{"chem"},$HMDB{$Data{$id}{"id"}},$Data{$id}{'Kingdom'},$Data{$id}{'Super_Class'},$Data{$id}{'Class'},$Data{$id}{'Sub_Class'},$Data{$id}{'Source'},$Data{$id}{'pathway'},$temp[1],$Data{$id}{"mz"};
                        }
                }
                if($count ==0)
                {
                        printf Out1 "%s\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%s\tNA\n",$temp[0],$temp[1];
                }
        }
}

close In1;
close Out1;

my %Gene;
open Ref4,"<",$dir0."/ALL_HMDB_gene.txt";
while(<Ref4>)
{
        chomp;
        my @temp=split /\s+/;
        $Gene{$temp[0]}{$temp[1]}=$temp[2];
}

open In2,"<",$dir."/".$fn.".HMDB.dat";
open Out2,">",$dir."/".$fn.".gene.dat";
while(<In2>)
{
        if(!/^ID/)
        {
                chomp;
                my @temp=split /\s+/;
                for my $n(sort {$a <=> $b} keys %{$Gene{$temp[2]}})
                {
                        printf Out2 "%s\t%s\t%s\t%s\t%s\n",$temp[0],$temp[1],$temp[2],$n,$Gene{$temp[2]}{$n};
                }
        }
}

close In2;
close Out2;

open In3,"<",$dir."/".$fn.".gene.dat";
open Out4,">",$dir."/".$fn.".gene.nred.1.dat";
open Out5,">",$dir."/".$fn.".gene.nred.2.dat";
my %Nred;
my %NRGene;
while(<In3>)
{
        chomp;
        my @temp=split /\s+/;
        $Nred{$temp[0]}{$temp[4]}=1;
        $NRGene{$temp[4]}=1;
}
for my $m (sort keys %Nred)
{
        my $count=1;
        for my $n (sort keys %{$Nred{$m}})
        {
		if(length($n)>0)
		{
                	printf Out4 "%s\t%d\t%s\n",$m,$count,$n;
                	$count++;
		}
        }
}
for my $i(sort keys %NRGene)
{
	if(length($i)>0)
	{
        	printf Out5 "%s\n",$i;
	}
}

close In3;
close Out4;
close Out5;
