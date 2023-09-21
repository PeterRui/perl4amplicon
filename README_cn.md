这套Perl脚本由兰州大学芮俊鹏博士编写。如果您有任何关于脚本的问题，请咨询芮博士 (peter_rjp@163.com)。
脚本需要使用命令行执行，例如Windows系统的命令提示符（CMD）、Linux系统的终端。
运行脚本时，如果不提供任何参数，会显示脚本的帮助信息（例如全部参数的说明）。

如有必要，请引用以下文章：
Rui J, Zhao Y, Cong N, Wang F, Li C, Liu X, Hu J, Ling N and Jing X (2023) Elevational distribution and seasonal dynamics of alpine soil prokaryotic communities. Front. Microbiol. 14:1280011. doi: 10.3389/fmicb.2023.1280011


1. perl trim_primer_in_fq.pl
     
该脚本用于检测fastq格式序列中的引物，在前后引物外侧的碱基（例如barcode）将被去除。建议去除引物，您也可以选择保留引物。

[多余的碱基][前引物][目的序列][后引物][多余的碱基]

用法: perl trim_primer_in_fq.pl -i [fastq文件] -l [引物列表文件] -d [输出目录] <可选参数>

-i的值可以是单个文件，也可以是若干文件，例如-i “fq/_*.fastq”表示fq目录下的所有扩展名为fastq的文件，使用通配符_*时请务必使用英文引号。

引物列表文件包含两行，依次是前引物、后引物序列及名称（请使用Tab制表符分隔），例如：

GMRCCIGGIGTIGGYTGYGC	nifH-2F

TTGTTGGCIGCRTASAKIGCCAT	nifH-3R


2. subsample_in_table.pl
     
该脚本用于对OTU/ASV表进行标准化，也就是说，将表格中的序列数目进行随机重抽样，使得每个样本中的序列总数保持一致。

注意: 输入文件每列为样本，每行为物种 (例如OTU、ASV)或功能，物种分类或其他额外信息可以放在输入文件的最后几列。

用法: perl subsample_in_table.pl -i [输入文件] -o [输出文件] -s [样本数目] <可选参数>


3. alpha_diversity.pl
     
该脚本基于物种丰度表（例如OTU表、ASV表）计算alpha多样性，例如Observed species, Chao1 richness, Shannon diversity, Simpson index, Pielou evenness, Good's coverage等。

注意: 输入文件每列为样本，每行为物种 (例如OTU、ASV)或功能，物种分类或其他额外信息可以放在输入文件的最后几列。

用法: perl alpha_diversity.pl -i [输入文件] -o [输出文件] -s [样本数目]


4. percent_in_table.pl
     
该脚本可将包含序列数目的物种丰度表（例如OTU表、ASV表）转化为百分数(不显示%)表格。

注意: 输入文件每列为样本，每行为物种 (例如OTU、ASV)或功能，物种分类或其他额外信息可以放在输入文件的最后几列。

用法: perl percent_in_table.pl -i [输入文件] -n [样本数目]


5. sum_taxa_from_otu_table.pl
     
该脚本基于物种丰度表（例如OTU表、ASV表）计算每个分类等级的物种总丰度，例如Clostridium属在每个样本中的丰度。
MG-RAST分级的的COG、Subsystems、KO功能表也可以作为输入文件，统计每个功能的丰度。 

用法: perl sum_taxa_from_otu_table.pl -i [输入文件] -s [样本数目] -D [数据类型] <可选参数>

-D: 输入文件的数据类型 (默认 0)，0代表百分比 (不含%，即percent_in_table.pl的输出文件)，1代表序列数目。

注意：物种分类请放在输入文件的最后几列（制表符分隔）。
例如界、门、纲、目、科、属六列，则对应的输出文件是taxa_L1.txt ~ taxa_L6.txt，它们不包含当前分类下的未分类序列丰度（因此每个样本的总和可能低于100%）；
如果您使用参数-a 1，会额外得到taxafull_L1.txt ~ taxafull_L6.txt，它们包含当前分类下的未分类序列丰度（也就是说，每个样本的总和为100%）。


6. rename_seqid_for_usearch.pl
     
该脚本用于fasta/fastq格式序列的重命名，ID重命名为"文件名_编号;barcodelabel=文件名;"。
来自不同样本的序列可合并到一个文件中进行分析（例如OTU聚类），Qiime和Usearch软件可识别序列的样本来源。

例如，文件"S1.fasta"的第一条序列重命名为"S1_1;barcodelabel=S1;"，第二条序列重命名为"S1_2;barcodelabel=S1;"。

用法: perl rename_seqid_for_usearch.pl -i [输入文件] -f [输入文件格式] <可选参数>

-i: 如有多个文件，可以使用通配符*（例如"*.fq"），或多个文件之间以英文分号隔开（例如"1.fq;2.fq;3.fq"），请务必使用英文引号！

-f: 输入文件的格式，0 - fasta，1 - fastq


7. fastq2fasta.pl
     
该脚本用于将fastq格式序列文件转为fasta格式文件。

Usage:perl fastq2fasta.pl [输入文件]
