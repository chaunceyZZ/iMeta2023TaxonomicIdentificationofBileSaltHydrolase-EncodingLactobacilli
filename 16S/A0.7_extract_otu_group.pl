use strict;
use File::Basename;

open Ref,"<","*.info";
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

#printf Out "ID\t";
#for my $i(sort keys %ID)
#{
#        printf Out "%s\t",$i;
#}
#printf Out "\n";

open In,"<",$InFile;
open Out,">",$OutFile;
while(<In>)
{
	if(/^ID/)
	{
		my @temp=split /\t/;
		if(exists $ID{$temp[0]})
		{               
		printf Out "%s\t%s\n",$temp[0];
		}	
	}
	
}
