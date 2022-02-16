import os

from Classes.DataForHTMLOutput import DataForHTMLOutput
from Classes.Gene import Gene
from Classes.Virus import Virus
from Classes.SamRecord import SamRecord
from Classes.Convertor import Convertor, create_html_output
from Classes.BlastApi import BlastApi

RNA_SAMPLE = "http://homel.vsb.cz/~vas218/files/viruses/viruses.filtered.sam"
VIRUS_NAME = "http://homel.vsb.cz/~vas218/files/viruses/viral_id.csv"
FEATURE_TO_GET_FROM_GTF = "cds"

if __name__ == '__main__':

    #BlastApi.send_query("blastn", "nt", "AGGATATTGTATTAGACCTGCAACCTCCAGACCCTGTAGGGTTACATTGCTATGAGCAATTAGTAGACAGCGCAGA")

    # Load gtf files
    genes, viruse_ids_with_no_gtf = Convertor.load_gtf_files_only_cds_or_gene(feature="CDS")



    # Creating map of VirusSeq -> VirusName
    map = Convertor.get_map_of_id_to_name(VIRUS_NAME)

    # Load a file and mapping ids to names, list of Sequences
    # [VirusSeq, Count] -> [VirusName, Count] | VirusSeq : VirusName = N : 1
    sam_records = Convertor.load_sam_files(RNA_SAMPLE, map, genes)

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
        Convertor.get_seqs_with_count_grouped_by(sam_records, "virus_name")

    count_of_mapped_records = len(mapped_records)
    count_of_not_mapped_records = len(not_mapped_records)

    # Take best candidates for Blast
    candidates = Convertor.choose_best_candidates_for_blast(not_mapped_records,
                                                            filter_repeating=False,
                                                            filter_acgt_probability=False,
                                                            filter_acgt_probability_from_fasta=True,
                                                            gtf_files=genes,
                                                            feature="CDS")
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
    genes_with_coverage = Gene.write_to_file_genes_with_percents(genes=genes, only_non_empty=True)

    # For each virus get % of mapped
    virus_with_coverage = Gene.write_to_file_virus_with_percents(genes=genes, only_non_empty=True)

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
                              sam_records_long_ends_starts= sam_records_long_ends_starts, map_virus_id_name=map)

    create_html_output(datas)