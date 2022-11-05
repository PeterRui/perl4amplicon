#!/usr/bin/perl -w

sub usage(){
<<EOF
This script is used to randomly subsample the number of reads in the OTU table.
Note: Each column is a sample. Each row is a(n) OTU/gene. It's OK to keep extra columns (e.g., taxonomy) in the OTU table.

Usage: perl subsample_in_table.pl -i [input file] -o [output file] -s [number of samples] <optional parameters> 

Parameters: 
-i the input file
-o the output file (default "otu_table_resample.txt")
-s number of samples
-n number of reads to be selected in each sample (optional, default using the number of the sample including the fewest reads) 
-d delimitor type in the OTU table (default 0): 0 - tab, 1 - comma
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

# === calculate the least reads ===
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
print "Sample $t[$min] contains the least reads ($sum[$min]). Resample to $n reads per sample\n";

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
