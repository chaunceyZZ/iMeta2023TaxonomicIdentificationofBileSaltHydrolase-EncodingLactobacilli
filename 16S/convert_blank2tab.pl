use strict;
my @files=glob("*.info");
for my $file(@files)
{
	open In,"<",$file;
	open Out,">","./tt/".$file;
	while(<In>)
	{
		chomp;
		my @temp=split  /\s+/;
		for my $i(0..$#temp)
		{
			printf Out "%s\t",$temp[$i];
		}
		printf Out "\n";
	}
}
