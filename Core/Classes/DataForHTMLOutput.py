from Core.Classes.Gene import Gene
from Core.Classes.SamRecord import SamRecord
from Core.Classes.Virus import Virus

#viruses_with_counts: list[list[Virus, count_of_all, count_of_amb, count_of_mapped_in, count_of_mapped_out]]

class DataForHTMLOutput:
    def __init__(self, viruses_with_counts,virus_coverage, genes_coverage, interesting_results_from_api,
                 total_count_of_sam_records, count_of_sam_records_with_no_gtf, count_of_sam_records_filtered_out_with_n,
                 count_of_mapped_in_sam_records, count_of_candidates_for_blast, count_of_filtered_out_candidates,
                 sam_records_long_ends_starts, map_virus_id_name, count_of_not_mapped_records,
                 find_sam_records_long_ends_starts, blast_is_set, best_candidates,
                 default_json_v, defaul_json_g, default_output
                 ): # type: (list[list[Virus, int, int, int, int]],list[list[int, float, int]],dict[str, list[Gene]], dict[str, int], int, int, int, int, int, int, list[SamRecord], dict[str: str], int, bool, bool, list[SamRecord], str, str, str) -> None

        self.virus_coverage = virus_coverage
        self.genes_coverage = genes_coverage
        self.viruses_with_counts = viruses_with_counts
        self.interesting_results_from_api = interesting_results_from_api
        self.total_count_of_sam_records = total_count_of_sam_records
        self.count_of_sam_records_with_no_gtf = count_of_sam_records_with_no_gtf
        self.count_of_sam_records_filtered_out_with_n = count_of_sam_records_filtered_out_with_n
        self.count_of_mapped_in_sam_records = count_of_mapped_in_sam_records
        self.count_of_not_mapped_records = count_of_not_mapped_records
        self.count_of_candidates_for_blast = count_of_candidates_for_blast
        self.count_of_filtered_out_candidates = count_of_filtered_out_candidates
        self.find_sam_records_long_ends_starts = find_sam_records_long_ends_starts
        self.sam_records_long_ends_starts = sam_records_long_ends_starts
        self.map_virus_id_name = map_virus_id_name
        self.blast_is_set = blast_is_set
        self.best_candidates = best_candidates
        self.default_json_v = default_json_v
        self.defaul_json_g = defaul_json_g
        self.default_output = default_output

