#!/usr/bin/perl -w

sub usage(){
<<EOF
This script can normalize the number of sequences in an ASV/OTU-table. In other words, the total numbers of reads in all samples will be the same after random subsample. The input file is an ASV/OTU-table with sequence numbers.

Note: Each column in the input table represents a sample, whereas each row represents a species (e.g., OTU/ASV). Additional columns in the input table, such as taxonomic infomation, are acceptable.

Usage: perl subsample_in_table.pl -i [input file] -o [output file] -s [number of samples] <optional parameters> 

Parameters: 
-i the input file
-o the output file (default "otu_table_resample.txt")
-s number of samples
-n number of reads to be kept in each sample (optional, using the sample size with the fewest reads as a default) 
-d delimitor type in the input table (default 0): 0 - tab, 1 - comma
-z how to treat the rows with no reads after subsample? (default 0)
	0 - remove them
	1 - keep them

EOF
}

use Getopt::Std;
use List::Util qw(shuffle);
getopts("i:o:s:n:d:z:",\%opts);

my $input= $opts{i};
my $out  = $opts{o}; $out = "otu_table_resample.txt" unless (defined($out));
my $s    = $opts{s};
my $n    = $opts{n}; 
my $delim= $opts{d};
my $zero = $opts{z};
print "Copyright: Junpeng Rui, Lanzhou University. peter_rjp\@163.com\n";
print "Please cite this article:\nRui J, Zhao Y, Cong N, Wang F, Li C, Liu X, Hu J, Ling N and Jing X (2023) Elevational distribution and seasonal dynamics of alpine soil prokaryotic communities. Front. Microbiol. 14:1280011. doi: 10.3389/fmicb.2023.1280011\n\n";
die usage() unless ($opts{i} and $opts{s});

print localtime()." - subsample in the otu table $input\n";

my ($d,$tit,$i,$j,$k,$tmp,$li,$min,$row,$sum_r,@f,@box,@sum,@a);

open TMP, $input || die "Cannot open the input file $input\n";
unless ($out) {
  $out=$input;
  $out =~ s/\.[\w]{1,3}$//;
  if ($delim) {
	$out .= "_resample.csv";
  } else {
	$out .= "_resample.txt";
  }
}
open OUT,">$out" || die "Cannot open the output file $out\n";

if ($delim) {
	$d=",";
} else {
	$d="\t";
}

# === read the OTU table ===
$tit=<TMP>;	# the title of OTU table
print OUT $tit;
$tit=~s/[\r\n]//g;
@t=split/$d/,$tit;
$i=0;
while ($li=<TMP>) {
	$li=~s/[\r\n]//g;
	@f=split/$d/,$li;
	$a[$i]=[@f];
	
	$j=1;
	while ($j<=$s) {
		$sum[$j]+=$f[$j];
		$j++;
	}
	$i++;
}
close TMP;
$row=$i;

# === look for the sample with the fewest reads ===
$min=1;
print "$t[1] contains $sum[1] reads.\n";
$j=2;
while ($j<=$s) {
	print "$t[$j] contains $sum[$j] reads.\n";
	if ($sum[$min]>$sum[$j]) {
		$min=$j;
	}
	$j++;
}
if (!$n or $sum[$min]<$n) {
	$n=$sum[$min];
}
print "Sample $t[$min] contains the fewest reads ($sum[$min]). Resample to $n reads per sample\n";

# === Sort Randomly ===
$j=1;
while ($j<=$s) {
	$i=0;
	@box=();
	while ($i<$row) {
		$k=0;
		while ($k<$a[$i][$j]) {
			push @box, $i;
			$k++;
		}
		$a[$i][$j]=0;
		$i++;
	}
	
	@box = shuffle(@box);	# sort randomly
	$k=0;
	while ($k<$n) {
		$a[$box[$k]][$j]++;
		$k++;
	}
	$j++;
}

# === Output the subsampled table ===
$i=0;
while ($i<$row) {
	if ($zero) {	# keep rows without reads
		$tmp=join $d,@{$a[$i]};
		print OUT "$tmp\n";
	} else {	# remove rows without reads
		$j=1;
		$sum_r=0;
		while ($j<=$s) {
			$sum_r+=$a[$i][$j];	# sum reads in the row i
			$j++;
		}
		if ($sum_r>0) {
			$tmp=join $d,@{$a[$i]};
			print OUT "$tmp\n";
		}
	}
	$i++;
}
close OUT;
