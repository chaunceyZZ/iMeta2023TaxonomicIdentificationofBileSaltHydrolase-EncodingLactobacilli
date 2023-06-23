use strict;
use File::Basename;
die("Argument: InFile RefFile OutFile\n") if ((@ARGV < 4) or (@ARGV >4));
my $Out=$ARGV[0];
my $file1=$ARGV[1];
my $file2=$ARGV[2];
my $flag=$ARGV[3];


my @files;
if($flag == 1){

	@files=glob($file1);
} else {
	my @file1=glob($file1);
	my @file2=glob($file2);
	@files=(@file1,@file2);
}
open Out1,">",$Out;
for my $file(@files)
{
	$file=~/(\w+)\.sub/;
	my $ID=$1;
	my $count=0;
	open In,"<",$file;

	while(<In>)
	{
		if(/^>/)
		{
			printf Out1 ">%s_%s\n",$ID,$count;
			$count++;
		}
		else
		{
			printf Out1 "%s",$_;
		}
	}
}
