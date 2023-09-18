#!/usr/bin/perl -w

sub usage(){
<<EOF
This script is used to rename sequence IDs in fasta/fastq-format files, and add ";barcodelabel=<File name>;" after the IDs, which could be identified by both USEARCH and Qiime. 
For example, the first ID in "S1.fasta" will named as "S1_1;barcodelabel=S1;", and the second ID will be "S1_2;barcodelabel=S1;".

Parameters: 
-i input files (e.g., "*.fasta" or "1.fq;2.fq;3.fq")
-s suffix of output (default is ".fasta" or ".fastq")
-d output dir (default "rename")
-f format of input files? (default 0)
	0 - fasta, 1 - fastq
EOF
}

use File::Basename qw<basename dirname>;
use File::Path;
use Getopt::Std;
getopts("i:s:d:f:",\%opts);

my $inputs= $opts{i}; 
my $suff  = $opts{s};
my $dout  = $opts{d}; $dout = "rename" unless (defined($dout));
my $format= $opts{f};
print "Copyright: Junpeng Rui, Lanzhou University. peter_rjp\@163.com\n";
print "Please cite this article:\nRui J, Zhao Y, Cong N, Wang F, Li C, Liu X, Hu J, Ling N and Jing X (2023) Elevational distribution and seasonal dynamics of alpine soil prokaryotic communities. Front. Microbiol. 14:1280011. doi: 10.3389/fmicb.2023.1280011\n\n";
die usage() unless ($opts{i});

unless (-e $dout) {
	mkpath($dout);
}

unless ($suff) {
  if ($format) {
	$suff=".fastq";
  } else {
	$suff=".fasta";
  }
}

my @myfile;
if ($inputs=~ /\*/) {  # wildcard mode
  @myfile = glob $inputs;
}
else {  # file list mode
  @myfile = split/;/,$inputs;
}
print "@myfile\n";

my $begin=1;

while(@myfile) {
	my $input = shift(@myfile);
	open MYIN, $input || die "Cannot open input file $input!\n";
	my $prefix=basename $input;
$prefix=~ s/\.(fq|fastq|fasta|fa|fas)$//i;
	open MYOUT, ">$dout/$prefix$suff" || die "Cannot open the output file.\n";
	print "$input -> $dout/$prefix$suff\n";
 
	my $n=0;
	$n+=$begin;
	my $nl=0;
	if ($format) {	# fastq
	  while($li=<MYIN>){
		$nl++;
		if ($nl % 4 == 1) {
			$li= '@'."${prefix}_$n;barcodelabel=${prefix};\n";
			$n++;
		}
		print MYOUT $li;
	  }
	} else {	# fasta
	  while($li=<MYIN>){
		if($li =~ /^\>/) {
			$li= '>'."${prefix}_$n;barcodelabel=${prefix};\n";
			$n++;
		}
		print MYOUT $li;
	  }
	}
	close MYOUT;
	close MYIN;
}


