from Classes.Gene import Gene
from Classes.SamRecord import SamRecord
from Classes.Virus import Virus

#viruses_with_counts: list[list[Virus, count_of_all, count_of_amb, count_of_mapped_in, count_of_mapped_out]]

class DataForHTMLOutput:
    def __init__(self, viruses_with_counts: list[list[Virus, int, int, int, int]],
                 virus_coverage: list[list[int, float, int]], genes_coverage: dict[str, list[Gene]], interesting_results_from_api: list[str],
                 total_count_of_sam_records: int, count_of_sam_records_with_no_gtf: int, count_of_sam_records_filtered_out_with_n: int,
                 count_of_mapped_in_sam_records: int, count_of_candidates_for_blast:int, count_of_filtered_out_candidates:int,
                 sam_records_long_ends_starts: list[SamRecord], map_virus_id_name: dict[str: str]):

        self.virus_coverage = virus_coverage
        self.genes_coverage = genes_coverage
        self.viruses_with_counts = viruses_with_counts
        self.interesting_results_from_api = interesting_results_from_api
        self.total_count_of_sam_records = total_count_of_sam_records
        self.count_of_sam_records_with_no_gtf = count_of_sam_records_with_no_gtf
        self.count_of_sam_records_filtered_out_with_n = count_of_sam_records_filtered_out_with_n
        self.count_of_mapped_in_sam_records = count_of_mapped_in_sam_records
        self.count_of_candidates_for_blast = count_of_candidates_for_blast
        self.count_of_filtered_out_candidates = count_of_filtered_out_candidates
        self.sam_records_long_ends_starts = sam_records_long_ends_starts
        self.map_virus_id_name = map_virus_id_name

