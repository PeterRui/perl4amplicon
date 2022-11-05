#!/usr/bin/perl -w
print "Copyright: Junpeng Rui, Lanzhou University. peter_rjp\@163.com\n";
print "This script is used to transfer fastq to fasta files.\n";
die "Usage:perl fastq2fasta.pl [the input file]\n" unless @ARGV>0;

my $fastq_file = $ARGV[0];
my $fasta_file = $fastq_file;
open TMP, "$fastq_file" || die "Can't open fastq $fastq_file!\n";    
$fasta_file =~ s/fastq$|fq$/fasta/; 
open MYOUT, ">$fasta_file"|| die "Can't open fasta $fasta_file"."_id!\n";
my ($li,$nl);
while($li=<TMP>){
$nl++;
if ($nl % 4 == 1) { # Illumina fastq contains 4 rows per seq. The qual data may also include @.
print MYOUT $li if $li =~ s/^\@/\>/; 
} 
elsif ($nl % 4 == 2) {
print MYOUT $li;  
}
}
close MYOUT;
close TMP;
print localtime()." - Completed: $fastq_file to $fasta_file\n";
