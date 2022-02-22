import argparse
import os
import random

from dominate.tags import a

from Core.Classes.DataForHTMLOutput import DataForHTMLOutput
from Core.Classes.Gene import Gene
from Core.Classes.Convertor import Convertor, create_html_output
from Core.Classes.BlastApi import BlastApi
from Core.Classes.Preprocess import preprocess

RNA_SAMPLE = "http://homel.vsb.cz/~vas218/files/viruses/viruses.filtered.sam"
VIRUS_NAME = "http://homel.vsb.cz/~vas218/files/viruses/viral_id.csv"

DEFAULT_VIRAL_ID = os.path.join("Data", "viral_id.csv")
DEFAULT_GTF_FOLDER = os.path.join("Data", "dir_gtf_files")
DEFAULT_FASTA_FOLDER = os.path.join("Data", "dir_fasta_genomes_files")
DEFAULT_REF_GENOME = os.path.join("Data", "ref_genome", "virus_genomes.fasta")
DEFAULT_JSON_VIRUSES = os.path.join("Output", "json", "virus_coverage_output.json")
DEFAULT_JSON_GENES = os.path.join("Output", "json", "genes_coverage_output.json")

def is_valid_file(filepath:str) -> bool:
    if not os.path.exists(filepath):
        return False
    else:
        return True


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="SAM file analyzer")
    parser.add_argument(action="store", dest="fastq_file", help="Input FASTQ file to align, filter and then analyze SAM file.", type=str)

    #parser.add_argument(action="store", dest="")

    """
    parser.add_argument("-v", "--viral_id_path", action="store", dest="viral_id_file",
                        help="Path to viral id file.", type=str, default=DEFAULT_VIRAL_ID)
    parser.add_argument("-gtf", "--gtf_folder", action="store", dest="gtf_folder",
                        help="Folder in which are folders with .gtf files. Look at the structure in Data/dir_gtf_files/",
                        type=str, default=DEFAULT_GTF_FOLDER)
    """
    parser.add_argument("-g", "--gtf_feature", action="store", dest="gtf_feature",
                        help="Take account only gtf records with this feature. 'CDS' or 'gene', default is 'CDS'.", default="CDS",
                        type=str)
    parser.add_argument("-r", "--recreate_ref_genome", action="store_true", dest="recreate_ref_genome",
                        help="If set, program will recreate and will create index for new referential genome from fasta files in Data/dir_fasta_genomes_files",
                        default=False)
    parser.add_argument("-e", "--empty", action="store_true", dest="empty",
                        help="If set, program will show in output.html all viruses and genes even if they have 0 %% coverage. Otherwise program will show only non 0 %% covered viruses and genes. ",
                        default=False)
    parser.add_argument("-L", "--find_longs", action="store_true", dest="find_longs",
                        help="Find records with long repeating substring on start or end of sequence", default=False)
    parser.add_argument("-l", "--find_longs_max_len", action="store", dest="l_n",
                        help="Set if argument -l is set. Specifies minimal length of substring for argument -l", default=25,
                        type=int)

    parser.add_argument("-B", "--blast", action="store_true", dest="blast", help="With this argument program will try to search blast api for some not mapped records. Warning! This can take even 5 minutes.",
                        default=False)
    parser.add_argument("-t", "--blast_time_limit", action="store", dest="blast_time_limit",
                        help="Time limit in seconds for search if Blast search is set.", default=5*60, type=int)
    parser.add_argument("-c", "--blast_max_candidates", action="store", dest="n_candidates",
                        help="Number of SAM records to be sent to Blast Api", default=10, type=int)
    parser.add_argument("-R", "--filter_repeating", action="store_true", dest="filter_repeating",
                        help="If set, program will filter SAM records for Blast api if atleast MIN_REPEATING %% of its sequence is repeating substring",
                        default=False)
    parser.add_argument("-x", "--filter_repeating_min_percent", action="store", dest="min_repeating", type=float, default=0.3,
                        help="Minimal [PERCENTS] %% of seq must be repeating substring to be filtered out. This argument can change [PERCENTS]. Enter percents in interval (0.0, 1.0), default value is 0.3")
    parser.add_argument("-P", "--filter_by_acgt_probability", action="store_true", dest="filter_acgt_probability",
                        help="If set, program will filter SAM records for Blast api, that doesn't have similar probabilities (0.25 +- deviation) of each nucleotide in sequence.",
                        default=False)
    parser.add_argument("-p", "--filter_by_acgt_probability_deviation", action="store", dest="dev_f2", type=float, default=0.05,
                        help="Deviation for -f2 filter option. Default value is 0.05")
    parser.add_argument("-O", "--filter_by_acgt_probability_from_fasta", action="store_true", dest="filter_acgt_probability_from_fasta",
                        help="If set, program will filter SAM records for Blast api, that doesn't have similar nucleotides probability with atleast one genome from fasta files",
                        default=False)
    parser.add_argument("-o", "--filter_by_acgt_probability_from_fasta_deviation", action="store", dest="dev_f3", default=0.05, type=float,
                        help="Deviation for -f3 filter option. Default value is 0.05")

    arguments = parser.parse_args()

    if not is_valid_file(arguments.fastq_file):
        print("The file %s does not exist!" % arguments.fastq_file)
        exit(-1)

    if arguments.l_n <= 0:
        print(f"Please enter argument -l1 as integer number greater than 0.")
        exit(-1)

    if arguments.dev_f3 <= 0:
        print(f"Please enter argument -f3D as float number greater than 0.")
        exit(-1)

    if arguments.dev_f2 <= 0:
        print(f"Please enter argument -f2D as float number greater than 0.")
        exit(-1)

    if arguments.n_candidates <= 0 or arguments.n_candidates > 100:
        print(f"Please enter argument -b2 as integer number greater than 0 and lower than or equal 100")
        exit(-1)

    if arguments.blast_time_limit <= 0:
        print(f"Please enter argument -b1 as integer number greater than 0.")
        exit(-1)

    if not (0.0 < arguments.min_repeating < 1.0):
        print(f"Please enter argument -f1M in interval (0.0, 1.0)")
        exit(-1)

    if arguments.gtf_feature not in ["CDS", "gene"]:
        print("GTF feature must be 'CDS' or 'gene'.")
        exit(-1)

    # Unzip, Align, Filter
    sam_file = preprocess(fastq_file=arguments.fastq_file, recreate_ref_genome=arguments.recreate_ref_genome,
                          DEFAULT_FASTA_FOLDER=DEFAULT_FASTA_FOLDER, DEFAULT_REF_GENOME=DEFAULT_REF_GENOME)

    # Load gtf files
    genes, viruse_ids_with_no_gtf = Convertor.load_gtf_files_only_cds_or_gene(directory=DEFAULT_GTF_FOLDER,
                                                                              feature=arguments.gtf_feature)

    # Creating map of VirusSeq -> VirusName
    map = Convertor.get_map_of_id_to_name(path=DEFAULT_VIRAL_ID)

    # Load a file and mapping ids to names, list of Sequences
    # [VirusSeq, Count] -> [VirusName, Count] | VirusSeq : VirusName = N : 1
    sam_records = Convertor.load_sam_files(sam_file, map, genes)

    total_count_of_sam_records = len(sam_records)

    # Get records with missing gtf file
    sam_records, sam_records_with_missing_gtf = Convertor.filter_records_with_missing_gtf(sam_records, viruse_ids_with_no_gtf)

    count_of_sam_records_with_no_gtf = len(sam_records_with_missing_gtf)

    # getting list [ Virus, count_of_all, count_of_amb ] by attribute of Sequence
    # Basic filtering
    count_of_sam_records_before_n_filter = len(sam_records)
    virus_with_count, sam_records_long_ends_starts, mapped_records, not_mapped_records =\
        Convertor.get_seqs_with_count_grouped_by(sam_records, "virus_name", feature=arguments.gtf_feature,
                                                 find_long_ends_starts=arguments.find_longs, max_len_of_end_start=arguments.l_n)

    count_of_mapped_records = len(mapped_records)
    count_of_not_mapped_records = len(not_mapped_records)

    # Work with BLAST API
    interesting_results_from_api = []
    best_candidates = []
    count_of_filtered_out_candidates = 0
    count_of_candidates = 0
    if arguments.blast:
        # Take best candidates for Blast
        candidates = Convertor.choose_best_candidates_for_blast(not_mapped_records,
                        filter_repeating=arguments.filter_repeating, max_percent_limit=arguments.min_repeating,
                        filter_acgt_probability=arguments.filter_acgt_probability, deviation_for_filtering_probability=arguments.dev_f2,
                        filter_acgt_probability_from_fasta=arguments.filter_acgt_probability_from_fasta, deviation_for_filtering_probability_from_fasta=arguments.dev_f3,
                        gtf_files=genes,
                        feature=arguments.gtf_feature,
                        directory_of_fastas=DEFAULT_FASTA_FOLDER)
        count_of_candidates = len(candidates)
        count_of_filtered_out_candidates = len(not_mapped_records) - len(candidates)

        if arguments.n_candidates > count_of_candidates:
            arguments.n_candidates = count_of_candidates
        best_candidates = random.sample(candidates, arguments.n_candidates)

        blast_results = BlastApi.send_multiple_queries("blastn", "nt", best_candidates, arguments.blast_time_limit)
        interesting_results_from_api = BlastApi.analyze_results_from_blast(blast_results)

    # For each gene get % of mapped
    genes_with_coverage = Gene.write_to_file_genes_with_percents(filename=DEFAULT_JSON_GENES, genes=genes,
                                                                 only_non_empty=not arguments.empty, map=map)

    # For each virus get % of mapped
    virus_with_coverage = Gene.write_to_file_virus_with_percents(filename=DEFAULT_JSON_VIRUSES, genes=genes,
                                                                 only_non_empty=not arguments.empty, map=map)

    # sort it by count desc
    virus_with_count.sort(key=lambda virus: virus[1], reverse=True)

    # Printing result VirusId\tVirusName\tCountOfAll\tCountOfAmb
    print("virus_id\tvirus_name\tcount_of_all\tcount_of_amb\tcount_of_IN\tcount_of_OUT\n")
    for item in virus_with_count:
        print(f"[SAM Analyzer]: {item[0].virus_id}\t{item[0].virus_name}\t{item[1]}\t{item[2]}\t{item[3]}\t{item[4]}")

    # Summary info
    print(f"\n[SAM Analyzer]: Total count of SAM records = {total_count_of_sam_records}")
    print(f"[SAM Analyzer]: Count of SAM records with missing GTF = {count_of_sam_records_with_no_gtf}")
    print(f"[SAM Analyzer]: Count of filtered out SAM records containing \"N\" = {count_of_sam_records_before_n_filter - (count_of_mapped_records + count_of_not_mapped_records)}")
    print(f"[SAM Analyzer]: Count of Mapped in SAM records = {count_of_mapped_records}")
    print(f"[SAM Analyzer]: Count of not Mapped in SAM records = {count_of_not_mapped_records}")
    if arguments.blast:
        print(f"[SAM Analyzer]: Count of Candidate SAM records for Blast api = {count_of_candidates}, only {len(best_candidates)} candidates were sent to api")
        print(f"[SAM Analyzer]: Count of filtered out Candidate SAM records for Blast api = {count_of_filtered_out_candidates}")
    if arguments.find_longs:
        print(f"[SAM Analyzer]: Count of SAM records that start or end with long repeating nucleotide [MIN_LENGTH={arguments.l_n}] = {len(sam_records_long_ends_starts)}")

    datas = DataForHTMLOutput(viruses_with_counts=virus_with_count, virus_coverage=virus_with_coverage, genes_coverage=genes_with_coverage,
                              blast_is_set=arguments.blast, best_candidates=best_candidates,
                              interesting_results_from_api=interesting_results_from_api, total_count_of_sam_records=total_count_of_sam_records,
                              count_of_sam_records_filtered_out_with_n=count_of_sam_records_before_n_filter - (count_of_mapped_records + count_of_not_mapped_records),
                              count_of_candidates_for_blast=count_of_candidates, count_of_filtered_out_candidates=count_of_filtered_out_candidates,
                              count_of_mapped_in_sam_records=count_of_mapped_records, count_of_sam_records_with_no_gtf=count_of_sam_records_with_no_gtf,
                              sam_records_long_ends_starts= sam_records_long_ends_starts, find_sam_records_long_ends_starts=arguments.find_longs,
                              map_virus_id_name=map, count_of_not_mapped_records=count_of_not_mapped_records)

    create_html_output(datas, sam_file)
