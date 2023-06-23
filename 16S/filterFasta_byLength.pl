#!/usr/bin/env perl
#
use File::Basename;
die("Argument: file.fa length\nKept seq >= length\n") if (@ARGV != 2);
my $fa = shift @ARGV;
my $len = shift @ARGV;
my($fan, $fadir, $suffix) = fileparse( $fa, ".fna");

my @output;
my ($seq_id, $seq);
open FA, "<$fa";
while(chomp($line = <FA>)){
	if($line =~ /^>/){
		push @output, "$seq_id\n$seq\n" if (length($seq) >= $len);
		$seq_id = $line;
		$seq = '';
	}else{
		$seq .= "$line";
	}
}
push @output, "$seq_id\n$seq\n" if (length($seq) > $len);
close FA;
my $output = $fadir . $fan . ".filter_${len}" . $suffix;
open LEN,">$output";
print LEN @output;
close LEN;
print "filter $fa by $len!\n";

