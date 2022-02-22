# sam-file-analyzer
#### Dependencies
- Python libraries: dominate (```pip install dominate```)
- Ubuntu tools: bwa (```apt install bwa```), samtools (```apt install samtools```), gunzip (```apt install gunzip```)
#### How to use + Examples
- In ubuntu run ```python3 main.py -h``` to see all possible arguments with description
- In ubuntu run ```python3 main.py inputFastQFile.fastq.gz``` for running tool with default parameters and without using Blast Api. Optimals Parameters are set.
- In ubuntu run ```python3 main.py inputFastQFile.fastq.gz -r``` if you want to recreate referential genome which will be in folder _Data/ref_genome/_ from all fasta files in folder _Data/dir_fasta_genomes_files/_
- In ubuntu run ```python3 main.py inputFastQFile.fastq.gz -B -c 20``` for running tool with default parameters and using Blast Api without filtering to search max 20 candidate sam files. 
- In ubuntu run ```python3 main.py inputFastQFile.fastq.gz -B -c 20 -R -P -O``` for running tool with using Blast api with all 3 filtering methodes with default parameters.
- #### Filtering candidates for Blast API
  - In ubuntu run ```python3 main.py inputFastQFile.fastq.gz -B -c 20 -R -x 0.3``` for filtering using only 1. filter methode (Filter repeating substring) with parameter -x 0.3 which means that SAM record that have atleast 30 % of their sequence covered by repeating substring will be filtred out of candidates for blast api. Parameter -x must be in range (0.0, 1.0)
  - In ubuntu run ```python3 main.py inputFastQFile.fastq.gz -B -c 20 -P -p 0.05``` for filtering using only 2. filter methode (Filter SAM records that doesn't have similar probability of each nucleotide in sequence in range 0.25 - deviation < PROB(A/C/G/T) < 0.25 + deviation). Parameter -p is the deviation and must be positive number.
  - In ubuntu run ```python3 main.py inputFastQFile.fastq.gz -B -c 20 -O -o 0.05``` for filtering using only 3. filter methode (Filter SAM records that doesn't have similar probability of nucleotides like atleast one genome in fasta files). Parameter -o is the deviation and must be positive number.
- In ubuntu run ```python3 main.py inputFastQFile.fastq.gz -e``` for running tool with default parameters and without Blast Api and with showing all viruses and genes even if they are not even partly covered.
- Main HTML output will be located in folder **Output/*.html**

#### Problems or questions 
- GTF files from NC_002016.1 to NC_002023.1 are same
- GTF files from NC_002204.1 to NC_002211.1 are same
- GTF files from NC_005218.1 to NC_005219.1 and NC_005222.1 are same
- GTF files from NC_006306.2 to NC_006312.2 are same
- GTF files from NC_007543.1 to NC_007547.1 and NC_007569.1 to NC_007574.1 are same
- GTF files from NC_011500.2 to NC_011510.2 are same
- GTF files from NC_021541_1 to NC_021551.1 are same
- !GTF file U21247.1 not found as assembly, only as Nucleotide
- !GTF files gi|9626032|lcl|HPV2REF.1|, gi|397005|lcl|HPV3REF.1|, gi|333074|lcl|HPV8REF.1|, gi|396910|lcl|HPV12REF.1|, gi|60295|lcl|HPV13REF.1|, gi|396918|lcl|HPV14REF.1|, gi|333245|lcl|HPV39REF.1| not found as assembly
- Problem when loading attributes which are delimited by ';' if semicolon is in value (...; note "Some words and; another words";... -> [..., 'note "Some words and', 'another words', ...]). Even when using csv.reader. Problem is that whole item isn't in " ", but only value ( key "value" -> to make it work -> "key value").

### Changelist:
* 27.10. 2021:
  * Runtime reduced from almost 2 mins to 25 secs by using dynamic programming when assigning genes to sam records.
  * Not loading whole gtf file, just gene and cds lines.
  * Added loading SEQ, CIGAR, POS from sam file to sam record.
  * Getting count_of_IN, count_of_OUT for CDS or GENE. But getting too many OUT, maybe error.
  * count getting by checking if POS from sam file is in <gene_start-78, gene_end) in atleast one gene associated with virus_id.
* 31.10. 2021:
  * Modified class Gene (added length_of_gene, records -> record, coverage_array)
  * Added calculating % of gene which are covered in sam file. (Gene.write_to_file_genes_with_percents)
  * Added calculating % of gene grouped by virus_id. (Gene.write_to_file_virus_with_percents)
  * ? Include or not the end index. Now not including.
* 08.11. 2021:
  * Modified Convertor.get_seqs_with_count_grouped_by() so it takes in account cigar string and takes mapped nucleotids only if cigar operation is "M". ???
  * Tried using Bio.Blast.NCBIWWW library with methode qblast(). Can return HTML, XML, Text, ASN.1. Accepts sequence ("AGGATATTGTATTAGACCTGCAACCTCCAGACCCTGTAGGGTTACATTGCTATGAGCAATTAGTAGACAGCGCAGA"), but it takes a very long time ( 5 - 10 minutes )
* 02.02. 2022:
  * Added filtering for all SAM records, which excludes records containing "N"
  * Added filtering for all SAM records, which excludes records without their gtf files
  * Added partition of SAM records to Mapped records and not mapped records
  * Virus and Genes coverage is calculating only for Mapped records.
  * Added option for getting list of SAM records (both mapped and not mapped) which contain A^n or C^n or G^n or T^n on start or end. Parameter n can change user.
  * Added 3 options to filter not mapped SAM records (Repeating substrings, Equal probability of Nucleotides (0.25 - deviation <= P(A,C,G,T) <= 0.25 + deviation), Checking if sam_record has similar Nucleotides probabilities like any virus from fasta files )
  * Zatím projdu celý soubor fasta a z něho celého počítám pravděpodobnosti. Možná lepší použít CDS/gene z gtf a počítat přímo z těch pozic????  
  * Added numerical summary (count of filtered, count of mapped, count of total, ...)
* 04.02. 2022:
  * Added option to use blast api for best 10 candidates. Using threads for making requests and getting results.  
  * Co s vystupem??
* 09.02. 2022:
  * Updated filtering option - "Checking if sam_record has similar Nucleotides probabilities like any virus from fasta files", checks only in areas (start - end) described by gtf files and the records CDS xor gene
  * U některých SAM souborů problém, že chybí gtf soubour, nebo chybí záznamy CDS xor gene... -> ignorování.
  * Updated work with api -> asking api about state every 60 seconds to comply with usage guidelines and avoid being blocked.
* 10.02. 2022:
  * Added time limit for queris on Blast api + After the time limit, program drops the query on api.
  * getting list of interesting results from results from blast api
* 16.02. 2022:
  * Modified result outputs from some functions
  * Main HTML output file was made. Contains Virus coverage table, Gene coverage table 
* 18.02. 2022:
  * Added graphs to html output using JS. Pie graphs show percentual coverage of each virus by its genes.
  * Made tables more user friendly (Show more, sorted records, sorted graphs)
  * .txt files with percentual coverage -> .json files, because of transfer of data to JS graphs.
* 22.02. 2022:
  * Tool is in final state
  * Tuned html output, Added preprocess part (Unzip, CreateRefGenome, Align, Filter), Added parameter parsing using argparse. 
