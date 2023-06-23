use strict;

use File::Basename;
die("Argument: InFile OutDir\n") if ((@ARGV < 2) or (@ARGV >2));

#my $InFile="low.LDA.Sig.mz.dat";
#my $dir="./";

my $InFile=$ARGV[0];
my ($fn, $dir1, undef) = fileparse($InFile, qw/.mz.dat/);
my $dir=$ARGV[1];



my %Data;
open Ref1,"<","/home/Metagenome-2/Ref_database/HMDB/megadb.20170119-121933.txt";
while(<Ref1>)
{
        my @temp=split /\t/;
        $Data{$temp[0]}{"mz"}=$temp[3];
        $Data{$temp[0]}{"chem"}=$temp[4];
}
open Ref2,"<","/home/Metagenome-2/Ref_database/HMDB/HMDB.20170119-121933.txt";
while(<Ref2>)
{
        my @temp=split /\t/;
        $Data{$temp[0]}{"id"}=$temp[1];
}
my %HMDB;
open Ref3,"<","/home/Metagenome-2/Ref_database/HMDB/ALL_HMDB_name.txt";
while(<Ref3>)
{
        my @temp=split /\t/;
        $HMDB{$temp[0]}=$temp[1];
}

open In1,"<",$InFile;
open Out1,">",$dir."/".$fn.".HMDB.dat";
printf Out1 "ID\tHit_num\tHMDB_ID\tChem\tClass\tObs_mz\tRef_mz\n";

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
                                printf Out1 "%s\t%d\t%s\t%s\t%s\t%s\t%s\n",$temp[0],$count,$Data{$id}{"id"},$Data{$id}{"chem"},$HMDB{$Data{$id}{"id"}},$temp[1],$Data{$id}{"mz"};
                        }
                }
                if($count ==0)
                {
                        printf Out1 "%s\t0\tNA\tNA\tNA\t%s\tNA\n",$temp[0],$temp[1];
                }
        }
}

close In1;
close Out1;

my %Gene;
open Ref4,"<","/home/Metagenome-2/Ref_database/HMDB/ALL_HMDB_gene.txt";
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
