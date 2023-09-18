#!/usr/bin/perl -w

sub Usage(){
<<EOF
Alpha diversity indices can be calculated with this script, such as Observed species, Chao1 richness, Shannon diversity, Simpson index, Pielou evenness, and Good's coverage. The input file is an ASV/OTU-table with sequence numbers.

Note: Each column in the input table represents a sample, whereas each row represents a species (e.g., OTU/ASV). Additional columns in the input table, such as taxonomic infomation, are acceptable.

Usage: perl alpha_diversity.pl -i otu_table_resample.txt -o alpha_div.txt -s 3

Parameters: 
-i the input file (default "otu_table_resample.txt")
-o the output file (default "alpha_div.txt")
-s number of samples

Example of the input file (otu_table_resample.txt):
SampleID	Sample1	Sample2	Sample3	taxonomy
OTU_1001	20111	23222	22315	Bacteria;Proteobacteria
OTU_1002	19802	22101	13418	Bacteria;Acidobacteria
...
Note: taxonomy is an optional extra column.

EOF
}

use Getopt::Std;
getopts("i:o:s:",\%opts);

my $input = $opts{i}; $input = "otu_table_resample.txt" unless (defined($input));
my $out   = $opts{o}; $out = "alpha_div.txt" unless (defined($out));
my $sam   = $opts{s};
print "Copyright: Junpeng Rui, Lanzhou University. peter_rjp\@163.com\n";
print "Please cite this article:\nRui J, Zhao Y, Cong N, Wang F, Li C, Liu X, Hu J, Ling N and Jing X (2023) Elevational distribution and seasonal dynamics of alpine soil prokaryotic communities. Front. Microbiol. 14:1280011. doi: 10.3389/fmicb.2023.1280011\n\n";
die Usage() unless $opts{s};

open TMP, $input || die "Cannot open the input file\n";
open OUT, ">$out" || die "Cannot open the output file\n";

my $li=<TMP>;
chomp $li;
my @f=split/\t/,$li;
my $title=join "\t", @f[1..$sam];
print OUT "\t$title\n";

my ($c,$row,$p,$ni,$N,$div,@sum,@a,@chao,@shan,@even,@simp1,@simp2,@o1,@o2,@otu,@cover);

my $r=1;	# row
while ($li=<TMP>) {
	next unless $li=~/\d/;
	chomp($li);
	@f=split/\t/,$li;
	$c=1;	# column
	while ($c<=$sam) {
		# for shannon and simpson:
		$sum[$c]+=$f[$c];	# sum of each column
		$a[$r][$c]=$f[$c];
		# for chao1 and Good's:
		if ($f[$c]>0) {
			$otu[$c]++;
			if ($f[$c]==1) {
				$o1[$c]++;	# singletons
			}
			elsif ($f[$c]==2) {
				$o2[$c]++;	# doubletons
			}
		}
		$c++;
	}
	$r++;
}
close TMP;
$row=$r-1;	# total rows of data
print "There are $row rows.\n";

$c=1;
while ($c<=$sam) {
	$r=1;
	my ($H,$D,$chao,$cover);
	$N =$sum[$c];
	while ($r<=$row) {
	  $ni=$a[$r][$c];
	  if ($ni and $N) {
		$p =  $ni/$N;
		$H += -$p*log($p);
		$D += $p*$p;
	  }
	  $r++;
	}
	$chao = $otu[$c]+$o1[$c]*($o1[$c]-1)/(2*($o2[$c]+1));
	$cover = (1-$o1[$c]/$N)*100;
	push @chao, sprintf("%d",$chao);
	push @shan, sprintf("%.3f",$H);
	push @even, sprintf("%.3f",$H/log($otu[$c]));
	push @simp1,sprintf("%.3f",1/$D);
	push @simp2,sprintf("%.4f",1-$D);
	push @cover,sprintf("%.1f",$cover);
	$c++;
}
my $totalseq=join "\t", @sum[1..$sam];
my $obs_spe =join "\t", @otu[1..$sam];
my $chao1   =join "\t", @chao;
my $shannon =join "\t", @shan;
my $evenness=join "\t", @even;
my $simpson1=join "\t", @simp1;
my $simpson2=join "\t", @simp2;
my $good_cov=join "\t", @cover;

$div ="Number of reads (N)\t$totalseq\n";
$div.="Observed species (S)\t$obs_spe\n";
$div.="Chao1 index\t$chao1\n";
$div.="Shannon diversity (H)\t$shannon\n";
$div.="Pielou evenness (E)\t$evenness\n";
$div.="Simpson index (1-D)\t$simpson2\n";
$div.="Reciprocal Simpson (1/D)\t$simpson1\n";
$div.="Good's coverage (%)\t$good_cov\n";
print OUT $div;
close OUT;
