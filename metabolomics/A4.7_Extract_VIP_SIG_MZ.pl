use strict;
use File::Basename;
die("Argument: RefFile InFile OutFile\n") if ((@ARGV < 3) or (@ARGV >3));
	
my $RefFile=$ARGV[0];
my $InFile=$ARGV[1];
my $OutFile=$ARGV[2];



my %data;

open Ref,"<",$RefFile;
while(<Ref>)
{
	my @temp=split /\s+/;
	if(!/^\ID/)
	{
		$data{$temp[0]}=1;
	}
}

open In,"<",$InFile;
open Out,">",$OutFile;
my $count=0;
while(<In>)
{
	if($count == 0)
	{
		printf Out "%s",$_;
	}
	else
	{
		my @temp=split /\s+/;
		if(exists $data{$temp[0]})
		{
			printf Out "%s",$_;
		}
	}
	$count++;
}
