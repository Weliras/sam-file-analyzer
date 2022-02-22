import argparse
import os



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
    parser.add_argument("-gtf_f", "--gtf_feature", action="store", dest="gtf_feature",
                        help="Take account only gtf records with this feature. 'CDS' or 'gene', default is 'CDS'.", default="CDS",
                        type=str)
    parser.add_argument("-r", "--not_recreate_ref_genome", action="store_true", dest="recreate_ref_genome",
                        help="If set, program will recreate and will create index for new referential genome from fasta files in Data/dir_fasta_genomes_files",
                        default=False)
    parser.add_argument("-l", "--find_longs", action="store_true", dest="find_longs",
                        help="Find records with long repeating substring on start or end of sequence", default=False)
    parser.add_argument("-l1", "--find_longs_max_len", action="store", dest="l_n",
                        help="Set if argument -l is set. Specifies minimal length of substring for argument -l", default=25,
                        type=int)
    parser.add_argument("-b", "--blast", action="store_true", dest="blast", help="With this argument program will try to search blast api for some not mapped records. Warning! This can take even 5 minutes.",
                        default=False)

    arguments = parser.parse_args()

    #if not is_valid_file(arguments.fastq_file):
    #    print("The file %s does not exist!" % arguments.sam_file)
    #    exit(-1)
    """
    if not is_valid_file(arguments.viral_id_file):
        print("The file %s does not exist!" % arguments.viral_id_file)
        exit(-1)
    if not is_valid_file(arguments.gtf_folder):
        print("The Folder %s does not exist!" % arguments.viral_id_file)
        exit(-1)
    """
    if arguments.gtf_feature not in ["CDS", "gene"]:
        print("GTF feature must be 'CDS' or 'gene'.")
        exit(-1)

    # Unzip, Align, Filter
    sam_file = preprocess(fastq_file=arguments.fastq_file, recreate_ref_genome=arguments.recreate_ref_genome,
                          DEFAULT_FASTA_FOLDER=DEFAULT_FASTA_FOLDER, DEFAULT_REF_GENOME=DEFAULT_REF_GENOME)

    # Load gtf files
    genes, viruse_ids_with_no_gtf = Convertor.load_gtf_files_only_cds_or_gene(directory=DEFAULT_GTF_FOLDER,
                                                                              feature=arguments.gtf_feature)


    exit(-222)

    # Creating map of VirusSeq -> VirusName
    map = Convertor.get_map_of_id_to_name(path=DEFAULT_VIRAL_ID)

    # Load a file and mapping ids to names, list of Sequences
    # [VirusSeq, Count] -> [VirusName, Count] | VirusSeq : VirusName = N : 1
    sam_records = Convertor.load_sam_files(sam_file, map, genes)

    total_count_of_sam_records = len(sam_records)

    sam_records, sam_records_with_missing_gtf = Convertor.filter_records_with_missing_gtf(sam_records, viruse_ids_with_no_gtf)

    count_of_sam_records_with_no_gtf = len(sam_records_with_missing_gtf)
    # for sam_record in sam_records:
    #    sam_record.virus.genes = [gene for gene in genes if gene.virus_id == sam_record.virus.virus_id]


    # check for same virus names but different virus ids
    # sam_records.append(SamRecord(Virus("NC_002023.1", "Influenza_A_virus")))
    # sam_records.append(SamRecord(Virus("NC_002021.1", "Influenza_A_virus")))

    # getting list [ Virus, count_of_all, count_of_amb ] by attribute of Sequence
    # Basic filtering
    count_of_sam_records_before_n_filter = len(sam_records)
    virus_with_count, sam_records_long_ends_starts, mapped_records, not_mapped_records =\
        Convertor.get_seqs_with_count_grouped_by(sam_records, "virus_name", feature=arguments.gtf_feature,
                                                 find_long_ends_starts=arguments.find_longs, max_len_of_end_start=arguments.l_n)

    count_of_mapped_records = len(mapped_records)
    count_of_not_mapped_records = len(not_mapped_records)

    # Take best candidates for Blast
    candidates = Convertor.choose_best_candidates_for_blast(not_mapped_records,
                                                            filter_repeating=False,
                                                            filter_acgt_probability=False,
                                                            filter_acgt_probability_from_fasta=True,
                                                            gtf_files=genes,
                                                            feature=arguments.gtf_feature)
    count_of_candidates = len(candidates)
    count_of_filtered_out_candidates = len(not_mapped_records) - len(candidates)

    best_candidates = candidates[:10]

    # Fill in 10 candidates
    """
    if len(best_candidates) < 10:
        l = len(best_candidates)
        not_mapped_without_best_candidates = not_mapped_records.copy()
        for c in best_candidates:
            not_mapped_without_best_candidates.remove(c)

        #best_candidates.append(not_mapped_without_best_candidates[:10-l])
        best_candidates += not_mapped_without_best_candidates[:10-l]
    """

    blast_results = []
    #blast_results = BlastApi.send_multiple_queries("blastn", "nt", best_candidates, 20 * 60)


    interesting_results_from_api = BlastApi.analyze_results_from_blast(blast_results)

    # For each gene get % of mapped
    genes_with_coverage = Gene.write_to_file_genes_with_percents(genes=genes, only_non_empty=True, map=map)

    # For each virus get % of mapped
    virus_with_coverage = Gene.write_to_file_virus_with_percents(genes=genes, only_non_empty=True, map=map)

    # sort it by count desc
    virus_with_count.sort(key=lambda virus: virus[1], reverse=True)

    # Printing result VirusId\tVirusName\tCountOfAll\tCountOfAmb
    print("virus_id\tvirus_name\tcount_of_all\tcount_of_amb\tcount_of_IN\tcount_of_OUT\n")
    for item in virus_with_count:
        print(f"{item[0].virus_id}\t{item[0].virus_name}\t{item[1]}\t{item[2]}\t{item[3]}\t{item[4]}")

    # Summary info
    print(f"\nTotal count of SAM records = {total_count_of_sam_records}")
    print(f"Count of SAM records with missing GTF = {count_of_sam_records_with_no_gtf}")
    print(f"Count of filtered out SAM records containing \"N\" = {count_of_sam_records_before_n_filter - (count_of_mapped_records + count_of_not_mapped_records)}")
    print(f"Count of Mapped in SAM records = {count_of_mapped_records}")
    print(f"Count of not Mapped in SAM records = {count_of_not_mapped_records}")
    print(f"Count of Candidate SAM records for Blast api = {count_of_candidates}")
    print(f"Count of filtrered out Candidate SAM records for Blast api = {count_of_filtered_out_candidates}")

    datas = DataForHTMLOutput(viruses_with_counts=virus_with_count, virus_coverage=virus_with_coverage, genes_coverage=genes_with_coverage,
                              interesting_results_from_api=interesting_results_from_api, total_count_of_sam_records=total_count_of_sam_records,
                              count_of_sam_records_filtered_out_with_n=count_of_sam_records_before_n_filter - (count_of_mapped_records + count_of_not_mapped_records),
                              count_of_candidates_for_blast=count_of_candidates, count_of_filtered_out_candidates=count_of_filtered_out_candidates,
                              count_of_mapped_in_sam_records=count_of_mapped_records, count_of_sam_records_with_no_gtf=count_of_sam_records_with_no_gtf,
                              sam_records_long_ends_starts= sam_records_long_ends_starts, map_virus_id_name=map, count_of_not_mapped_records=count_of_not_mapped_records)

    create_html_output(datas)