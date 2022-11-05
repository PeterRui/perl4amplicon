#!/usr/bin/perl

sub Usage(){
<<EOF
This script is used to convert sequence-number table (e.g., an OTU table) to relative abundance table. 

Parameters: 
-i input files (e.g., "otu*.txt" or otu_table.txt)
-n number of samples (default value = column numbers - 1)
	Note: you should use -n if the table includes extra columns such as taxonomy.
-T total reads in each sample (optional)
-s suffix of output files (default "ab")
	e.g, The name of output files will be "***_ab.txt" if using -s "ab"
-z remove a row if the sum value of the rows is zero (default 0):
	0 - no, 1 - yes
EOF
}

use File::Basename;
use Getopt::Std;
getopts("i:n:T:s:z:",\%opts);

my $inputs= $opts{i};
my $numS  = $opts{n}; 
my $numT  = $opts{T};
my $suffix= $opts{s}; $suffix="ab" unless(defined($suffix));
my $rm0   = $opts{z};
print "Copyright: Junpeng Rui, Lanzhou University. peter_rjp\@163.com\n";
die Usage() unless ($opts{i});

my @sumT;
my @myfile;
if ($inputs=~ /\*/) {  # wildcard mode
	@myfile = glob $inputs;
} else {  # single file mode
	push @myfile,$inputs;
}
print "The input files are: @myfile\n";

my $sout="_";
if ($suffix) {
	$sout.=$suffix;
}

my $nf=0;
while (@myfile) {
  my ($input,$fname,$li,$n,$i,$j,$t1,$t2,$t,$dout,@a,@c,@o);  
  $input = shift @myfile;
  $fname = basename($input);
  $dout  = dirname($input);
  $fname =~ s/\.tsv$|\.txt$//i;

  open TMP, $input || die;
  open OUT, ">$dout/$fname$sout.txt" || die;

# print the title
  $li=<TMP>;
  print OUT $li;
# read the matrix
  $i=1;
  while($li=<TMP>) {
	chomp $li;
	@c=split /\t/,$li;
	$n=@c;
	if($i==1 and !$numS and !$nf) {
		$numS=$n-1;
	}
	$a[$i][0]=$c[0];  # otu id
	
	$j=1;
	while($j<=$numS) {
		$a[0][$j]+=$c[$j];  # sum reads of a sample
		$a[$i][$j]=$c[$j];  # store otu numbers into a(i,j)
		$j++;
	}
	while($j<$n) {
		$o[$i].="\t$c[$j]";  # other info such as taxa
		$j++;
	}
	$i++;
  }
  $otu=$i-1;
  close TMP;
  
  $j=1;
  while($j<=$numS) {
	$a[0][$j]=$numT if ($numT);
	$j++;
  }
  
# print the percent matrix
  $i=1;
  while($i<=$otu) {
	my $sumR;
	$j=1;
	$t=$a[$i][0];

	while($j<=$numS) {
	  $sumR+=$a[$i][$j];
	  if($a[0][$j]) {
		if($a[$i][$j]) {
			$t1=$a[$i][$j]/$a[0][$j]*100;
			$t2=sprintf "%.6f",$t1;
		}
		else {
			$t2=0;
		}
		$t.="\t$t2";
		$j++;
	  }
	  else {	# sum of a column is 0.
		$t.="\t0";
		$j++;
	  }
	}
	unless ($rm0 and $sumR==0) {
		print OUT "$t$o[$i]\n";
	}
	$i++;
  }
  close OUT;
  $nf++;
}
