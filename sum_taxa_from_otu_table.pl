#!/usr/bin/perl -w

sub Usage(){
<<EOF
Based on an ASV/OTU-table, this script can calculate the total abundance of each taxon at each taxonomic rank.
It is also appropriate for hierarchical tables like COG, subsystems, and KO in MG-RAST. 
Note: Values from different unclassified taxa will not be added up.

Parameters: 
-i input files
-o the output directory (default "taxa")
-s number of samples
-t number of taxa rank (default 6)
-D data type (default 0)
	0 - percentage (%); 1 - number of reads
-p the prefix of the output file (default "taxa")

Example of the OTU table (named "otu_table_resample_ab.txt"):
OTU_ID	Sam1	Sam2	Sam3	Domain	Phylum	Class	Order	Family	Genus
OTU_1	0.4523	1.3044	0.8532	Bacteria	Verrucomicrobia	Spartobacteria	Chthoniobacterales
OTU_2	0.5686	0.7810	0.7662	Bacteria	Gemmatimonadetes	Gemmatimonadetes	Gemmatimonadales	Gemmatimonadaceae	Gemmatimonas
OTU_3	0.3097	0.2650	0.1810	Bacteria	Proteobacteria	Betaproteobacteria	Burkholderiales	Alcaligenaceae
...

Example of usage: perl sum_taxa_from_otu_table.pl -i otu_table_resample_ab.txt -o out -s 3 -t 6 -D 0 -p taxa
the output files with unclassified taxa: out/taxafull_L{NUM}.txt
the output files without unclassified taxa: out/taxa_L{NUM}.txt
NUM=1~6 (-t 6, from Domain to Genus)

EOF
}

use Getopt::Std;
use File::Path;
getopts("i:o:s:t:D:p:",\%opts);

my $inputs= $opts{i};
my $dout  = $opts{o}; $dout="taxa" unless (defined($dout));
my $numS  = $opts{s};
my $numT  = $opts{t}; $numT = 6 unless (defined($numT));
my $dtype = $opts{D};
my $pout  = $opts{p}; $pout="taxa" unless (defined($pout));
print "Copyright: Junpeng Rui, Lanzhou University. peter_rjp\@163.com\n";
print "Please cite this article:\nRui J, Zhao Y, Cong N, Wang F, Li C, Liu X, Hu J, Ling N and Jing X (2023) Elevational distribution and seasonal dynamics of alpine soil prokaryotic communities. Front. Microbiol. 14:1280011. doi: 10.3389/fmicb.2023.1280011\n\n";
die Usage() unless ($opts{i} and $opts{s});

my @myfile;
if ($inputs=~ /\*/) {  # wildcard mode
	@myfile = glob $inputs;
}
else {  # single file mode
	push @myfile,$inputs;
}
print "@myfile\n";

my ($lv,$fname,@rk,@c);
mkpath($dout);

while (@myfile) {
  $input = shift @myfile;
  @rk=(1..$numT);
  
  while (@rk) {
	my ($li,$la,$i,$j,$k,$x,$y,$otu,$ln,$row,@a,@b,%raw,%ltaxa);
	$lv=shift @rk;
	$y=$numS+$lv;
	
	my $out1="$dout/${pout}full_L$lv.txt";
	my $out2="$dout/${pout}_L$lv.txt";
	open TMP, $input || die "Cannot open [$input].\n";
	open OUT, ">$out1" || die "Cannot open [$out1].\n";	# full data
	open OUT2, ">$out2" || die "Cannot open [$out2].\n";	# remove unclassified data

	# print the title
	$li=<TMP>;
	$li=~s/[\r\n]//g;
	@c=split /\t/,$li;
	$ln="Taxa";
	$j=1;
	while($j<=$y) {
		$ln.="\t$c[$j]";
		$j++;
	}
	print OUT "$ln\n";
	print OUT2 "$ln\n";

	# read the matrix
	$i=1;
	while($li=<TMP>) {
		my ($id,$lt);
		$li=~s/[\r\n]//g;
		@c=split /\t/,$li;
		$j=1;
		while($j<=$numS) {
			$a[$i][$j]=$c[$j];  # store data into a(i,j)
			$j++;
		}
		while($j<$y) {
			$id.="$c[$j]\t";	# taxa
			$lt=$c[$j] if $c[$j];	# the rank is not empty
			$j++;
		}
		$id.=$c[$j];
		$a[$i][0]=$id;
		$raw{$a[$i][0]}=$id;	# store the raw taxa
		
		if ($c[$j]) {
			$ltaxa{$a[$i][0]}=$c[$j];	# store the last rank of taxa
		} else {	# the last rank is empty
			$ltaxa{$a[$i][0]}="unclassified ".$lt;
		}
		
		$i++;
	}
	$otu=$i-1;
	
	# deal with the same taxa
	$j=1;
	while($j<=$numS) {
		$i=1;
		my %u;
		while($i<=$otu) {
			$u{$a[$i][0]}+=$a[$i][$j];  # sum of the same taxa
			$i++;
		}
		
		$x=1;
		foreach $k (sort keys %u) {
			if($u{$k}==0) {
				$b[$x][$j]=0;
			}
			else {
				if ($dtype) {
					$b[$x][$j]=$u{$k};
				} else {
					$b[$x][$j]=sprintf "%.4g",$u{$k};  # store sum data into b(x,j)
				}
			}
			$b[$x][0]=$k;
			$x++;
		}
		$row=$x-1;
		
		$j++;
	}
	print "$input\t$otu\t$row\n";
	
	# output the results
	$i=1;
	while($i<=$row) {
		$j=1;
		$ln=$ltaxa{$b[$i][0]};
		
		while($j<=$numS) {
			$ln.="\t".$b[$i][$j];
			$j++;
		}
		$ln.="\t".$raw{$b[$i][0]};
		print OUT "$ln\n";
		if (!($ltaxa{$b[$i][0]}=~/unclassified /)) {
			print OUT2 "$ln\n";
		}
		$i++;
	}
	close OUT;
	close OUT2;
	close TMP;
  }
}
