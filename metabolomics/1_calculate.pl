use strict;
my %Hits;
my @files=glob("YH.data.txt");
for my $file(@files)
{
    my $count; 
	my @temp;
    open In,"<",$file;
    open Out,">","1.".$file;
	while(<In>)
	{
		chomp;
		if(/^name/)
		{
			printf Out "%s\n", $_;
		}
		@temp=split /\s+/;
		$count=0;
		for (my $i=1;$i<@temp;$i++)
		{
			if($temp[$i]==0)
			{
				$count=$count+1;
			}
		}
		if($count<(@temp/2))
		{
			printf Out "%s\n", $_;
            $Hits{$temp[0]}=1;
		}
        if($count>(@temp/2))
        {
            #printf "%s\n", $temp[0];
        }
	}
}
open In1,"<","YH.mz.dat";
open Out1,">","1.YH.mz.dat";
print Out1 "name\tmzmed\tRtmed\n";
while(<In1>)
{
    chomp;
    my @temp = split /\s+/;
    if(exists $Hits{$temp[0]}){
        printf Out1 "%s\n", $_;
    }
}
