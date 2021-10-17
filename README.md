# sam-file-analyzer
TO-DO list
- Loading SAM file and reading virus_id - Done
- Mapping virus_id -> virus_name - Done
- Getting output in format Virus, count - Done
- Getting output in format Virus, count_of_all, count_of_amb - Done
- Download GTF files - 90 % Done
- Load GTF files to memory (new class, load files, check) - Done
- Look at gene position, count of inside/outside - Not Done

Problems or questions 
- GTF files from NC_002016.1 to NC_002023.1 are same
- GTF files from NC_002204.1 to NC_002211.1 are same
- GTF files from NC_005218.1 to NC_005219.1 and NC_005222.1 are same
- GTF files from NC_006306.2 to NC_006312.2 are same
- GTF files from NC_007543.1 to NC_007547.1 and NC_007569.1 to NC_007574.1 are same
- GTF files from NC_011500.2 to NC_011510.2 are same
- GTF files from NC_021541_1 to NC_021551.1 are same
- !GTF file U21247.1 not found as assembly, only as Nucleotide
- !GTF files gi|9626032|lcl|HPV2REF.1|, gi|397005|lcl|HPV3REF.1|, gi|333074|lcl|HPV8REF.1|, gi|396910|lcl|HPV12REF.1|, gi|60295|lcl|HPV13REF.1|, gi|396918|lcl|HPV14REF.1|, gi|333245|lcl|HPV39REF.1| not found as assembly
- Runtime with assigning genes to viruses - 1m 10s
