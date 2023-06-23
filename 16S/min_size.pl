#!/usr/bin/env perl
#
use File::Basename;
die("Argument: file.fa length\nKept seq >= length\n") if (@ARGV != 2);
my $fa = shift @ARGV;
my $len = shift @ARGV;




open FA, "<$fa";
while(chomp($line = <FA>)){
	if($line =~ /^>/){
		$seq = '';
	}else{
		$seq .= "$line"; ## 合并字符串
		if (length($seq) <= $len){
			$len=length($seq);
			
		}
	}
}
close FA;

print "min size $len\n";
