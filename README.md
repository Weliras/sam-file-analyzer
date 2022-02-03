# sam-file-analyzer
#### TO-DO list
- Loading SAM file and reading virus_id - Done
- Mapping virus_id -> virus_name - Done
- Getting output in format Virus, count - Done
- Getting output in format Virus, count_of_all, count_of_amb - Done
- Download GTF files - 90 % Done
- Load GTF files to memory (new class, load files, check) - Done
- Look at gene position, count of inside/outside - Done
- Get percents of coverage - Done
- Make code readible - Not done

#### Current state
- Got in memory all sam records with info about mapped virus and other ambiquous viruses. Every virus has associated genes (by virus_id).
- Getting output in format Virus, count_of_all, count_of_amb, count_of_IN, count_of_OUT + Percentages of coverege in files

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

#### Bugs
- If any attribute while reading GTF has semicolon in value -> Value is wrong. 

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
