# perl4amplicon
These Perl scripts are designed by Dr. Junpeng Rui. If you have any questions about these scripts, please email him (peter_rjp@163.com).
The scripts should be executed using command-line tools such as CMD on Windows or Terminal on Linux. 
Execute the scripts without any parameters will yield help information (e.g., optional parameters).

Please cite this article:
Rui J, Zhao Y, Cong N, Wang F, Li C, Liu X, Hu J, Ling N and Jing X (2023) Elevational distribution and seasonal dynamics of alpine soil prokaryotic communities. Front. Microbiol. 14:1280011. doi: 10.3389/fmicb.2023.1280011


  1. perl trim_primer_in_fq.pl
     
This script is used to search primers in fastq files. Bases outside them (e.g., barcodes) and/or primers will be trimmed:
[Bases to be trimmed][primer 1][good sequence][primer 2][Bases to be trimmed]

Usage: perl trim_primer_in_fq.pl -i [fastq files] -l [the primer list] -d [the output directory] <optional parameters>

If you use -i “fq/*.fastq”, it means all fastq files in the input directory “fq”.

Example of the text file including 2 primers:
GMRCCIGGIGTIGGYTGYGC	nifH-2F
TTGTTGGCIGCRTASAKIGCCAT	nifH-3R


  2. subsample_in_table.pl
     
This script can normalize the number of sequences in an ASV/OTU-table. In other words, the total numbers of reads in all samples will be the same after random subsample.
The input file is an ASV/OTU-table with sequence numbers.
Note: Each column in the input table represents a sample, whereas each row represents a species (e.g., OTU/ASV). Additional columns in the input table, such as taxonomic infomation, are acceptable.

Usage: perl subsample_in_table.pl -i [input file] -o [output file] -s [number of samples] <optional parameters>


  3. alpha_diversity.pl
     
Alpha diversity indices can be calculated with this script, such as Observed species, Chao1 richness, Shannon diversity, Simpson index, Pielou evenness, and Good's coverage. The input file is an ASV/OTU-table with sequence numbers.
Note: Each column in the input table represents a sample, whereas each row represents a species (e.g., OTU/ASV). Additional columns in the input table, such as taxonomic infomation, are acceptable.

Usage: perl alpha_diversity.pl -i [input file] -o [output file] -s [number of samples]


  4. percent_in_table.pl
     
This script will convert sequence numbers in an ASV/OTU table into percentages (%).
Note: Each column in the input table represents a sample, whereas each row represents a species (e.g., OTU/ASV). Additional columns in the input table, such as taxonomic infomation, are acceptable.

Usage: perl percent_in_table.pl -i [input file] -n [number of samples]


  5. sum_taxa_from_otu_table.pl
     
Based on an ASV/OTU-table, this script can calculate the total abundance of each taxon at each taxonomic rank.
It is also appropriate for hierarchical tables like COG, subsystems, and KO in MG-RAST. 
Note: Values from different unclassified taxa will not be added up.

Usage: perl sum_taxa_from_otu_table.pl -i [input file] -s [number of samples] -D [data type] <optional parameters>

-D: data type of the input file (default 0). 0 means relative abundance (%, i.e., percentage, the output file of percent_in_table.pl), while 1 means sequence numbers.


  6. rename_seqid_for_usearch.pl
     
This script is used to rename sequence IDs in fasta/fastq-format files, and add ";barcodelabel=<File name>;" after the IDs, which could be identified by both USEARCH and Qiime. 
For example, the first ID in "S1.fasta" will named as "S1_1;barcodelabel=S1;", and the second ID will be "S1_2;barcodelabel=S1;".

Usage: perl rename_seqid_for_usearch.pl -i [input files] -f [format of input files] <optional parameters>

-i: e.g., "*.fasta" or "1.fq;2.fq;3.fq"
-f: 0 - fasta, 1 - fastq


  7. fastq2fasta.pl
     
This script is used to convert fastq-format reads to fasta-format reads.

Usage:perl fastq2fasta.pl [input file]


