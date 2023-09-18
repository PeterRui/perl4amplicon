#!/usr/bin/perl -w

sub Usage(){
<<EOF
This script is used to search primers in fastq files. Bases outside them (e.g., barcodes) and/or primers will be trimmed:
[Bases to be trimmed][primer 1][good sequence][primer 2][Bases to be trimmed]

Usage: perl trim_primer_in_fq.pl -i [fastq files] -l [the primer list] -d [the output directory] <optional parameters>

Parameters: 
-i the input files in fastq format (e.g., -i "raw/*.fastq", or -i raw/1.fq)
-l a text file including primers
-d the output directory, which should be different from the input directory.
-r remove primers? (default 1)
	0 - No, only remove the bases outside them;
	1 - Yes, remove both the primers and the bases outside them.
-s keep reads which include one primer? (default 0)
	0 - No;
	1 - Yes.
-c remove sequences with more than one primer 1 or primer 2, such as primer-dimers? (default 1)
	0 - No;
	1 - Yes.

Example of the tab-delimited text file including 2 rows (i.e., 2 primers):
GMRCCIGGIGTIGGYTGYGC	nifH-2F
TTGTTGGCIGCRTASAKIGCCAT	nifH-3R

Note:The 1st column is primer sequence. The rest columns are useless. Degenerate bases are OK.
EOF
}

use Getopt::Std;
use File::Basename;
use File::Path;
getopts("i:l:d:r:s:c:",\%opts);

my $inputs= $opts{i};
my $list  = $opts{l};
my $dout  = $opts{d};
my $rmPrm = $opts{r}; $rmPrm=1 unless defined($rmPrm);
my $onePrm= $opts{s};
my $cmr   = $opts{c}; $cmr=1 unless defined($cmr);
print "Copyright: Junpeng Rui, Lanzhou University. peter_rjp\@163.com\n";
print "Please cite this article:\nRui J, Zhao Y, Cong N, Wang F, Li C, Liu X, Hu J, Ling N and Jing X (2023) Elevational distribution and seasonal dynamics of alpine soil prokaryotic communities. Front. Microbiol. 14:1280011. doi: 10.3389/fmicb.2023.1280011\n\n";
die Usage() unless ($opts{i} and $opts{l} and $opts{d});
mkpath($dout);
my ($li,$id,$rc,@f,@cp,@rp,@lp);

my @myfile;
if ($inputs=~ /\*/) {  # wildcard mode
  @myfile = glob $inputs;
} else {  # single file mode
  push @myfile,$inputs;
}
print "@myfile\n";

# ===== read the text file including primers =====
open TMP, $list || die "Cannot open the text file [$list]\n";
while ($li=<TMP>) {
	$li=~s/[\r\n]//g;
	@f=split/\t/,$li;
	# check the primer
	my $tmp=$f[0];
	$tmp=~s/[ATGCRYMKBDHVSWNI]//gi;
	die "[Error] Unrecognized characters are detected: \"$tmp\"\n" if length($tmp);
	# store the primers
	if ($f[0]) {
		push @cp,$f[0];
		push @lp,length($f[0]);	# length of the primer
		# primer reverse complement
		$rc=reverse($f[0]);
		$rc=~tr/ATGCRYMKBDHVatgcrymkbdhv/TACGYRKMVHDBtacgtrkmvhdb/;
		push @rp,$rc;
	}
}
print "primer 1\t$cp[0]\t$rp[0]\n";
print "primer 2\t$cp[1]\t$rp[1]\n";

die "[Error] Require 2 primers!\n" if @cp<2;
warn "[Warning] More than 2 primers! The top 2 will be used!\n" if @cp>2;
close TMP;

my $log="trim.primer.log";
open LOG, ">$log" || die "Cannot open the log file $log\n";
$tmp="Files\tNum.of.reads\tWith.2.primers\tOnly.1st.primer\tOnly.2nd.primer\n";
print LOG $tmp;
print $tmp;

# ===== read the input files ======
while (@myfile) {
my $input=shift @myfile;
my $din=dirname($input);
die "[Error] The output directory should be different from the input directory.\n" if ($din eq $dout);
my $fname=basename($input);
my $out="$dout/$fname";

open TMP,$input || die "Cannot open the fastq file $input\n";
open OUT,">$out" || die "Cannot open the output file $out\n";

my $n_cmr=$n=$nR=$nD=$n1=$n2=0;
my ($redun,%num,%seq,%qlt);

while($li=<TMP>) {
	if ($n%4==0 && $li=~/^@(.*)\s*$/) {
		$id=$1;
		$num{$id}++;
		if ($num{$id}>1){
			$nR++;
			$redun.="$id\t$num{$id}\n";
		}
	} elsif ($n%4==1) {
		$li=~s/[\r\n]//g;
		#$li=~s/-|~|\s|\.//g;
		$seq{$id}=$li;
	} elsif ($n%4==3) {
		$li=~s/[\r\n]//g;
		$qlt{$id}=$li;
	}
	$n++;
}
close TMP;

foreach $id (sort keys %seq) {
	my ($new,$newq,$qrc);
	my $query=$seq{$id};
	my $qual =$qlt{$id};
	
	if($cmr) {	# check redundancy primers
		my $numP1=&number_of_primer($query,$cp[0])+&number_of_primer($query,$rp[0]);
		my $numP2=&number_of_primer($query,$cp[1])+&number_of_primer($query,$rp[1]);
		if ($numP1>1 or $numP2>1) {
			$n_cmr++;
			next;
		}
	}
	my ($f0,$f1)=&search_primer($query,$cp[0]);	# search 1st primer
	my ($r0,$r1)=&search_primer($query,$rp[1]);	# search 2nd primer (reverse complement)
	if ($f0==-1 and $r0==-1) {
		$qrc=reverse($qual);
		$rc=reverse($query);
		$rc=~tr/ATGCRYMKBDHVatgcrymkbdhv/TACGYRKMVHDBtacgtrkmvhdb/;	# reverse complement of the query sequence
		($f0,$f1)=&search_primer($rc,$cp[0]);
		($r0,$r1)=&search_primer($rc,$rp[1]);
		next if ($f0==-1 and $r0==-1);	# skip if failed
		$query=$rc;	# use reverse complement sequence if succeed
		$qual=$qrc;
	}
	if ($rmPrm) {	# remove primers
		if ($f0>=0 and $r0>$f1) {	# with 2 primers
			$new=substr($query,$f1+1,$r0-$f1-1);
			$newq=substr($qual,$f1+1,$r0-$f1-1);
			$nD++;
		} elsif ($onePrm) {	# also keep reads with only one primer
			if ($f0>=0 and $r0==-1) {
				$new=substr($query,$f1+1);
				$newq=substr($qual,$f1+1);
				$n1++;
			} elsif ($r0>0 and $f0==-1) {
				$new=substr($query,0,$r0);
				$newq=substr($qual,0,$r0);
				$n2++;
			}
		}
	} else {	# keep primers
		if ($f0>=0 and $r0>$f1) {	# with 2 primers
			$new=substr($query,$f0,$r1-$f0+1);
			$newq=substr($qual,$f0,$r1-$f0+1);
			$nD++;
		} elsif ($onePrm) {	# also keep reads with only one primer
			if ($f0>=0 and $r0==-1) {
				$new=substr($query,$f0);
				$newq=substr($qual,$f0);
				$n1++;
			} elsif ($r0>0 and $f0==-1) {
				$new=substr($query,0,$r1+1);
				$newq=substr($qual,0,$r1+1);
				$n2++;
			}
		}
	}
	
	print OUT "@"."$id\n$new\n+\n$newq\n" if ($new);
}
close OUT;
print "[Warning] $n_cmr reads are removed from [$out] due to more than one primer 1 or primer 2!\n" if $n_cmr;

my $nS=$n/4;
$tmp="$input\t$nS\t$nD\t$n1\t$n2\n";
print LOG $tmp;
print $tmp;

if ($nR>0) {
	my $fr=$input;
	$fr=~s/(fq|fastq)$//i;
	$fr.="redun.txt";
	open RD,">$fr" || die "Cannot open $fr\n";
	print RD "Redundant ID\tNum\n$redun";
	close RD;
	print "$nR IDs are redundant in [$input]! See the file [$fr]\n";
}

}
close LOG;

# ========= functions =========
sub IUB_to_regexp {
# Usage: Translate IUB ambiguity codes to regular expressions
my ($iub) = @_;
my $re = '';
my %iub2chr = (
	A => 'A',
	C => 'C',
	G => 'G',
	T => 'T',
	R => '(A|G)',
	Y => '(C|T)',
	M => '(A|C)',
	K => '(G|T)',
	S => '(C|G)',
	W => '(A|T)',
	B => '(C|G|T|S|K|Y)',
	D => '(A|G|T|R|K|W)',
	H => '(A|C|T|M|Y|W)',
	V => '(A|C|G|M|S|R)',
	N => '(A|C|G|T|R|Y|M|K|S|W|B|D|H|V|N|I)',
	I => '(A|C|G|T|R|Y|M|K|S|W|B|D|H|V|N|I)'
);
$iub =~ s/\^//g;
for ( my $i = 0 ; $i < length($iub) ; ++$i ) {
	$re.= $iub2chr{substr($iub, $i, 1)};
}
return $re;
}

sub search_primer {
# Usage: Search a primer in a sequence, and return the location.
my ($seq,$primer) = @_;
	my $exp=&IUB_to_regexp($primer);
	if ($seq=~/($exp)/i) {
		$c0=index($seq,$1);
		$c1=$c0+length($1)-1;
		return($c0,$c1);
	} else {
		return(-1,-1);	# failed
	}
}

sub number_of_primer {
# Usage: Search a primer in a sequence, and return the hit number.
my ($seq,$primer) = @_;
	my $exp=&IUB_to_regexp($primer);
	my $hit=($seq=~s/($exp)//gi);
	return($hit);
}