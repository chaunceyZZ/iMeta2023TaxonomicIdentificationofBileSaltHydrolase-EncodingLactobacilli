#!/usr/bin/env perl
#
use File::Basename;
die("Argument: _otu.txt _rep_set.fa cut-off[<4]\n") if ((@ARGV < 2) or (@ARGV >3));
my $otuFile = shift @ARGV;
my $repSet = shift @ARGV;
my $cut_off = 4;
$cut_off = (shift @ARGV) if (@ARGV);
my ($fn1, $dir1, undef) = fileparse($otuFile, qw/.txt/);
my ($fn, $dir, undef) = fileparse($repSet, qw/.fa/);

my %otu2reads_num;
my @output1;
open OTU,"<$otuFile";
while($line = <OTU>){
	my @line = split /\s+/, $line;
	next if (@line < ($cut_off + 1)); 
	if(!exists $otu2reads_num{$line[0]}){
		$otu2reads_num{$line[0]} = @line - 1;
		push @output1, $line;
	}else{
		die("Duplicated OTU id in $otuFile\n");
	}
}
close OTU;

my @output;
my $flag = 0;
open REPSET,"<$repSet";
while($line = <REPSET>){
	if($line =~ /^>/){
		my $otuid = (split /\s+/, $line)[0];
		$otuid =~ s/^>//;
		if(exists $otu2reads_num{$otuid}){
			$flag = 1;
			push @output, $line;
		}else{
			$flag = 0;
		}
	}else{
		push @output, $line if $flag;
	}
}
close REPSET;

my $output1 = $dir1 . $fn1 . ".filter_" . $cut_off . ".txt";
open OUTPUT, ">$output1";
print OUTPUT @output1;
close OUTPUT;

my $output = $dir . $fn . ".filter_" . $cut_off . ".fa";
open OUTPUT,">$output";
print OUTPUT @output;
close OUTPUT;
print "Filtered!\n";

