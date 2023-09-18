#!/usr/bin/perl -w
print "Copyright: Junpeng Rui, Lanzhou University. peter_rjp\@163.com\n";
print "Please cite this article:\nRui J, Zhao Y, Cong N, Wang F, Li C, Liu X, Hu J, Ling N and Jing X (2023) Elevational distribution and seasonal dynamics of alpine soil prokaryotic communities. Front. Microbiol. 14:1280011. doi: 10.3389/fmicb.2023.1280011\n\n";
print "This script is used to convert fastq-format reads to fasta-format reads.\n";

die "Usage:perl fastq2fasta.pl [the input file]\n" unless @ARGV>0;

my $fastq_file = $ARGV[0];
my $fasta_file = $fastq_file;
open TMP, "$fastq_file" || die "Cannot open the fastq file: $fastq_file!\n";    
$fasta_file =~ s/fastq$|fq$/fasta/i; 
open MYOUT, ">$fasta_file"|| die "Cannot open the fasta file!\n";
my ($li,$nl);

while($li=<TMP>){
	$nl++;
	if ($nl % 4 == 1) { # Illumina fastq contains 4 rows per seq. 
		$li =~ s/^@/\>/;
		print MYOUT $li;
	} elsif ($nl % 4 == 2) {
		print MYOUT $li;  
	}
}
close MYOUT;
close TMP;
print localtime()." - Completed: $fastq_file to $fasta_file\n";
