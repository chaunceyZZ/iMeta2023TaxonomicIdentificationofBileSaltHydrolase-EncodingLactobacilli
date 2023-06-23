use strict;
use File::Basename;
die("Argument: InFile grpFile OutFile groups\n") if ((@ARGV < 4) );

my $InFile=$ARGV[0];
my $grpFile=$ARGV[1];
my $OutFile=$ARGV[2];

my %Groups;
for my $i(3..$#ARGV)
{
	$Groups{$ARGV[$i]}=1;
}

open Ref,"<",$grpFile;
my %ID;
while(<Ref>)
{
	chomp;
	my @temp=split /\t/;
	if(exists $Groups{$temp[1]})
	{
		$ID{$temp[0]}=1;
	}
}
open Out,">",$OutFile;
printf Out "Phyla Genera\t";
for my $i(sort keys %ID)
{
        printf Out "%s\t",$i;
}
printf Out "\n";

open In,"<",$InFile;
my @key;
while(<In>)
{
	if(/^Phyla Genera/)
	{
		my @temp=split /\t/;
		for my $i(1..$#temp)
		{
			if(exists $ID{$temp[$i]})
			{
				@key=(@key,$i);
			}
		}
	}else{
		my @temp=split /\t/;
		printf Out "%s\t",$temp[0];
		for my $i(@key)
                {
                        printf Out "%s\t",$temp[$i];
                }
		printf Out "\n";
	}
	
}
