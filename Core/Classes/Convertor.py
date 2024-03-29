import datetime
import os
import shutil
import subprocess
import sys
import traceback
from copy import deepcopy
import dominate

from Core.Classes.DataForHTMLOutput import DataForHTMLOutput
from Core.Classes.Gene import Gene, GTF_File_Line
from Core.Classes.Virus import Virus
from Core.Classes.SamRecord import SamRecord
from Core.Classes.FastaFile import FastaFile
from dominate.tags import *


class Convertor:

    @staticmethod
    def get_map_of_id_to_name(path): # type: (str) -> dict[str, str]
        virusSeq_to_name_map = {}
        try:
            with open(path, "r") as file:
                for line in file:
                    #decoded_line = line.decode("utf-8")
                    decoded_line = line.strip()
                    line_split = decoded_line.split('\t')
                    virusSeq_to_name_map[line_split[0]] = line_split[1]
                    # print(decoded_line)
                virusSeq_to_name_map['*'] = '*'  # For unknown
        except Exception as e:
            print(e)
            traceback.print_exc(file=sys.stdout)
            exit(-1)
        return virusSeq_to_name_map

    @staticmethod
    def load_fasta_files(gtf_files, directory_of_fastas, feature = "CDS"): # type: (list[Gene], str, str) -> list[FastaFile]

        # GTF File preformat
        virus_id_to_gtf = dict()
        for gtf_file in gtf_files:
            if gtf_file.virus_id not in virus_id_to_gtf.keys():
                virus_id_to_gtf[gtf_file.virus_id] = []
            virus_id_to_gtf[gtf_file.virus_id].append(gtf_file)

        fasta_files = []
        directory = directory_of_fastas

        folders = [os.path.join(directory, o) for o in os.listdir(directory) if
                   os.path.isdir(os.path.join(directory, o))]
        for dir in folders:
            virus_name = dir.split("/")[2]
            fasta_file = [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f))]
            if len(fasta_file) <= 0:
                print(f"Fasta file for virus_name: {virus_name} was not found.")
                continue
            fasta_file = fasta_file[0]

            virus_id = str()
            sequence = str()
            nucleotides_probability_fasta = {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0}
            nucleotides_counts_fasta = {"A": 0, "C": 0, "G": 0, "T": 0}
            total_nucleotides_count = 0

            # reading fasta file, getting virus_id, sequence.
            try:
                with open(os.path.join(dir, fasta_file), "r") as file:
                    for line in file:
                        if line.startswith(">"):
                            virus_id = line.split()[0][1:].replace("|", "-")
                        else:
                            stripped_line = line.strip()
                            # total_nucleotides_count += len(stripped_line)
                            sequence += stripped_line
            except Exception as e:
                print(e)
                traceback.print_exc(file=sys.stdout)
                exit(-1)

            # Check for viruses that doesnt have gtf files or have but dont have records with feature CDS xor GENE
            if virus_id not in virus_id_to_gtf.keys():
                continue

            # Computing count of all nucleotides in seq by gtf end, start
            total_nucleotides_count_in_range_by_gtf = 0
            for gene in virus_id_to_gtf[virus_id]:
                if gene.record.feature == feature:
                    total_nucleotides_count_in_range_by_gtf += len(sequence[gene.record.start:gene.record.end])

            # Computing count of each nucleotide for all genes and then calculating probability
            for nucleotide in nucleotides_probability_fasta.keys():
                for gene in virus_id_to_gtf[virus_id]:
                    if gene.record.feature == feature:
                        nucleotides_counts_fasta[nucleotide] += sequence[gene.record.start:gene.record.end].count(
                            nucleotide)
                        # nucleotides_probability_fasta[nucleotide] = sequence[gene.record.start:gene.record.end].count(nucleotide)
                        # nucleotides_probability_fasta[nucleotide] = sequence.count(nucleotide) / total_nucleotides_count
                nucleotides_probability_fasta[nucleotide] = nucleotides_counts_fasta[nucleotide] / (
                    1 if total_nucleotides_count_in_range_by_gtf <= 0 else total_nucleotides_count_in_range_by_gtf)
            fasta_files.append(FastaFile(virus_id, virus_name, sequence, nucleotides_probability_fasta))
        return fasta_files

    @staticmethod
    def load_gtf_files_only_cds_or_gene(directory, feature="CDS"):  # type: (str, str) -> (list[Gene], list[str])
        '''
        Function returns only cds and gene lines from gtf files
        :return: list of genes and list of virus ids witch no gtf files
        '''
        d = directory
        folders = [os.path.join(d, o) for o in os.listdir(d) if os.path.isdir(os.path.join(d, o))]
        genes = []
        viruse_ids_with_no_gtf = []
        for dir in folders:
            #print(dir)
            #virus_id = dir.split("\\")[1].replace("-", "|")
            virus_id = dir.split("/")[2].replace("-", "|")
            gtf_file = [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f))]
            if len(gtf_file) <= 0:
                print(f"[SAM Analyzer]: Gtf file for virus_id: {virus_id} was not found.")
                viruse_ids_with_no_gtf.append(virus_id)
                continue
            gtf_file = gtf_file[0]
            try:
                with open(os.path.join(dir, gtf_file), "r") as file:
                    for line in file:
                        if line.startswith("#"):
                            continue

                        gtf_line = GTF_File_Line()
                        line_split = line.split("\t")
                        gtf_line.load_from_line(line_split)

                        if gtf_line.feature != feature:
                            continue
                        # if gtf_line.feature != "gene" and gtf_line.feature != "CDS":
                        #    continue
                        gene = Gene(virus_id)
                        gene.records.append(deepcopy(gtf_line))  # deprecated
                        gene.record = deepcopy(gtf_line)
                        genes.append(deepcopy(gene))

                        del gene
            except Exception as e:
                print(f"[SAM Analyzer]: {e}")
                traceback.print_exc(file=sys.stdout)
                exit(-1)
        return genes, viruse_ids_with_no_gtf

    @staticmethod
    def load_gtf_files(gtf_folder):  # type: (str) -> list[Gene]
        d = gtf_folder
        folders = [os.path.join(d, o) for o in os.listdir(d) if os.path.isdir(os.path.join(d, o))]

        genes = []
        for dir in folders:
            virus_id = dir.split("\\")[1].replace("-", "|")
            gtf_file = [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f))]
            if len(gtf_file) <= 0:
                print(f"[SAM Analyzer]: Gtf file for virus_id: {virus_id} was not found.")
                continue
            gtf_file = gtf_file[0]
            try:
                with open(os.path.join(dir, gtf_file), "r") as file:
                    prev_protein_id = ""
                    prev_gene_id = ""
                    first_in_file = True
                    gene = Gene(virus_id)
                    for line in file:

                        if line == "###\n":  # End of file, insert gene and delete
                            genes.append(deepcopy(gene))
                            del gene
                        if line.startswith("#"):
                            continue

                        gtf_line = GTF_File_Line()
                        line_split = line.split("\t")
                        gtf_line.load_from_line(line_split)

                        # if feature is transcript or exon it just appends them, because these features doesnt have ids
                        if gtf_line.feature == "transcript" or gtf_line.feature == "exon":
                            gene.records.append(gtf_line)
                            first_in_file = False
                            continue

                        if gtf_line.attributes["gene_id"] != "":
                            if not first_in_file:
                                if prev_gene_id != gtf_line.attributes["gene_id"]:
                                    genes.append(deepcopy(gene))
                                    del gene
                                    gene = Gene(virus_id=virus_id)
                        else:
                            if not first_in_file:
                                if prev_protein_id != gtf_line.attributes["protein_id"]:
                                    genes.append(deepcopy(gene))
                                    del gene
                                    gene = Gene(virus_id=virus_id)

                        # getattr(gene, gtf_line.feature).append(gtf_line)
                        gene.records.append(gtf_line)
                        first_in_file = False

                        if "gene_id" in gtf_line.attributes.keys():
                            prev_gene_id = gtf_line.attributes["gene_id"]
                        if "protein_id" in gtf_line.attributes.keys():
                            prev_protein_id = gtf_line.attributes["protein_id"]

            except Exception as e:
                print(f"[SAM Analyzer]: {e}")
                traceback.print_exc(file=sys.stdout)
                exit(-1)
        return genes

    @staticmethod
    def load_sam_files(path, map, genes):  # type: (str, dict[str, str], list[Gene]) -> list[SamRecord]
        """
        Using dynamic programming to get runtime under 30 secs
        :param genes: Loaded genes from GTF files
        :param path: Path to input SAM file
        :param map: Map of Virus id -> Virus Name
        :return: list of SamRecords loaded from a file. SamRecord->(Virus->Genes,AmbViruses->Genes)
        """
        try:
            dynamic_programming = dict()
            for virus_id in map.keys():
                dynamic_programming[virus_id] = None

            with open(path, "r") as file:
                sam_records = []
                for line in file:
                    #decoded_line = line.decode("utf-8")
                    line_split = line.split('\t')

                    ambiguous_viruses = []
                    for column in line_split:
                        if column.startswith("XA:Z:"):
                            ambiguous_viruses_column = column.split(":")[2].split(";")
                            for ambiguous_virus in ambiguous_viruses_column:
                                if ambiguous_virus != "" and ambiguous_virus != "\n":
                                    v = ambiguous_virus.split(",")
                                    if dynamic_programming[v[0]] is not None:
                                        ambiguous_viruses.append(Virus(virus_id=v[0], virus_name=map[v[0]],
                                                                       genes=dynamic_programming[v[0]]))
                                        continue
                                    # genes to virus
                                    g = [g for g in genes if g.virus_id == v[0]]
                                    ambiguous_viruses.append(Virus(virus_id=v[0], virus_name=map[v[0]], genes=g))
                                    dynamic_programming[v[0]] = g

                    if dynamic_programming[line_split[2]] is not None:
                        sam_records.append(
                            SamRecord(virus=Virus(virus_id=line_split[2], virus_name=map[line_split[2]],
                                                  genes=dynamic_programming[line_split[2]]),
                                      ambiguous_viruses=ambiguous_viruses, pos=line_split[3], cigar=line_split[5],
                                      seq=line_split[9], qname=line_split[0], flag=line_split[1], rname=line_split[2],
                                      mapq=line_split[4]))
                        continue
                    g = [g for g in genes if g.virus_id == line_split[2]]
                    sam_records.append(
                        SamRecord(virus=Virus(virus_id=line_split[2], virus_name=map[line_split[2]], genes=g),
                                  ambiguous_viruses=ambiguous_viruses, pos=line_split[3], cigar=line_split[5],
                                  seq=line_split[9], qname=line_split[0], flag=line_split[1], rname=line_split[2],
                                  mapq=line_split[4]))
                    dynamic_programming[line_split[2]] = g
        except Exception as e:
            print(f"[SAM Analyzer]: {e}")
            traceback.print_exc(file=sys.stdout)
            exit(-1)
        return sam_records

    @staticmethod
    def filter_records_with_missing_gtf(sam_records, viruse_ids_with_no_gtf):  # type: (list[SamRecord], list[str]) -> (list[SamRecord], list[SamRecord])
        missing_gtf = [rec for rec in sam_records if rec.virus.virus_id in viruse_ids_with_no_gtf]
        not_missing_gtf = [rec for rec in sam_records if rec.virus.virus_id not in viruse_ids_with_no_gtf]
        return not_missing_gtf, missing_gtf

    @staticmethod
    def get_seqs_with_count_grouped_by(sam_records, attr, feature="CDS",find_long_ends_starts=False, max_len_of_end_start=10):  # type: (list[SamRecord], str, str, bool, int) -> (list[Virus, int, int, int, int], list[SamRecord])
        """
        :param feature: "CDS" or "gene"
        :param sam_records: Loaded lines by method load_file()
        :param attr: Attribute of Virus class which is used to group viruses
        :param find_long_ends_starts: Find seqs with ends or starts containing A^n, C^n, G^n, T^n
        :param max_len_of_end_start: Length of found A^n or C^n or G^n or T^n on start or end of seq
        :return: (list[Virus, count_of_all, count_of_amb, count_of_mapped_id, count_of_mapped_out], grouped by attr), list of records which has A^n or C^n or G^n or T^n on start or end of seq
        """
        viruses_with_count = []
        sam_records_long_ends_starts = []
        count_of_all = len(sam_records)
        count_of_not_mapped = 0
        count_of_mapped = 0
        count_of_filtered_out = 0
        count_of_not_filtered_out = 0
        sam_records_mapped_in = []
        sam_records_mapped_out = []
        for sam_record in sam_records:

            # Filtering based on seq containing N + complex filtering
            if not filter_seq(sam_record, False):
                count_of_filtered_out += 1
                continue
            # Filtering based on repeating nucleotides
            # if filter_repeating and not calc_longest_repeating_subsequence(sam_record.SEQ, max_percent_limit, max_repeating_size):
            #    count_of_filtered_out += 1
            #    continue

            count_of_not_filtered_out += 1

            # Founding seq with T^n / A^n / C^n / G^n on start or end.
            if find_long_ends_starts:
                calc_long_ends_starts(sam_record, sam_records_long_ends_starts, max_len_of_end_start)

            count_of_genes_in_which_is_record = 0
            count_of_genes_in_which_is_not_record = 0

            # Decoding CIGAR string and computing count of nucleotides which were mapped in and out.
            for gene in sam_record.virus.genes:
                if gene.record.feature == feature:
                    count_of_genes_in_which_is_record, count_of_genes_in_which_is_not_record = \
                        calc_coverage_array_for_gene(count_of_genes_in_which_is_record,
                                                     count_of_genes_in_which_is_not_record, gene=gene,
                                                     sam_record=sam_record)

            # New Virus found and grouping
            if not any(s for s in viruses_with_count if getattr(s[0], attr) == getattr(sam_record.virus, attr)):
                tmp = [sam_record.virus, 1, 0, 0, 0]
                if count_of_genes_in_which_is_record > 0:
                    tmp[3] = 1
                elif count_of_genes_in_which_is_not_record > 0:
                    tmp[4] = 1
                viruses_with_count.append(tmp)
            # Virus with this ID already added
            else:
                ind = [viruses_with_count.index(item) for item in viruses_with_count if getattr(item[0], attr) ==
                       getattr(sam_record.virus, attr)]
                if count_of_genes_in_which_is_record > 0:
                    viruses_with_count[ind[0]][3] += 1  # Add count of IN
                elif count_of_genes_in_which_is_not_record > 0:
                    viruses_with_count[ind[0]][4] += 1  # Add count of OUT
                viruses_with_count[ind[0]][1] += 1

            # counting occurences of virus in other amb viruses
            count_of_amb_genes_mapped_in = 0
            for amb_virus in sam_record.ambiguous_viruses:
                count_of_amb_genes_in_which_is_record = 0
                count_of_amb_genes_in_which_is_not_record = 0

                # Calculating coverage array of amb viruses.
                for gene in amb_virus.genes:
                    if gene.record.feature == feature:

                        count_of_amb_genes_in_which_is_record, count_of_amb_genes_in_which_is_not_record = \
                            calc_coverage_array_for_gene(count_of_amb_genes_in_which_is_record
                                                         , count_of_amb_genes_in_which_is_not_record,
                                                         gene=gene, sam_record=sam_record)

                if not any(s for s in viruses_with_count if getattr(s[0], attr) == getattr(amb_virus, attr)):
                    tmp = [amb_virus, 1, 1, 0, 0]
                    if count_of_amb_genes_in_which_is_record > 0:
                        tmp[3] = 1
                    elif count_of_amb_genes_in_which_is_not_record > 0:
                        tmp[4] = 1
                    viruses_with_count.append(tmp)
                else:
                    ind = [viruses_with_count.index(item) for item in viruses_with_count if getattr(item[0], attr) ==
                           getattr(amb_virus, attr)]
                    viruses_with_count[ind[0]][1] += 1  # count of all
                    viruses_with_count[ind[0]][2] += 1  # count of amb
                    if count_of_amb_genes_in_which_is_record > 0:
                        viruses_with_count[ind[0]][3] += 1  # count of amb IN
                    elif count_of_amb_genes_in_which_is_not_record > 0:
                        viruses_with_count[ind[0]][4] += 1  # count of amb OUT
                if count_of_amb_genes_in_which_is_record > 0:
                    count_of_amb_genes_mapped_in += 1

            # If virus and other amb viruses haven't been mapped in, they are not mappped
            if count_of_amb_genes_mapped_in == 0 and count_of_genes_in_which_is_record == 0:
                count_of_not_mapped += 1
                sam_record.mapped = False
                sam_records_mapped_out.append(sam_record)
            else:
                count_of_mapped += 1
                sam_record.mapped = True
                sam_records_mapped_in.append(sam_record)



        return viruses_with_count, sam_records_long_ends_starts, sam_records_mapped_in, sam_records_mapped_out

    @staticmethod
    def choose_best_candidates_for_blast(sam_records_not_mapped, directory_of_fastas, filter_repeating=False, max_percent_limit=0.30, max_repeating_size=3,
                                         filter_acgt_probability=False, deviation_for_filtering_probability=0.05, filter_acgt_probability_from_fasta=False,
                                         gtf_files=None, feature="CDS", deviation_for_filtering_probability_from_fasta=0.05):  # type: (list[SamRecord], str, bool, float, int, bool, float, bool, list[Gene], str, float) -> list[SamRecord]
        """
        :param filter_repeating: Filter repeating seqs or not. If true -> calculations will take more time
        :param max_percent_limit: If filter_repeating is True, this percent means if n % is covering the seq and is greater then this max limit, then is this record considered invalid and is filtered out.
        :param max_repeating_size: Length of repeating Substring
        :return: list of candidates
        """

        # Load fasta files
        fasta_files = None
        if filter_acgt_probability_from_fasta:
            fasta_files = Convertor.load_fasta_files(gtf_files, directory_of_fastas, feature)

        count_of_filtered_out = 0
        candidates = []
        for sam_record in sam_records_not_mapped:

            # Filtering based on seq containing N
            if not filter_seq_with_n(sam_record):
                count_of_filtered_out += 1
                continue

            # Filtering based on ACGT probability, 0.25 - deviation <= P(A) == P(C) == P(G) == P(T) <= 0.25 + deviation
            if filter_acgt_probability and not filter_seq_with_different_probability_of_nucleotides(sam_record,
                                                                                                    deviation_for_filtering_probability):
                count_of_filtered_out += 1
                continue

            # Filtering based on ACGT probability from similar fasta files.
            if filter_acgt_probability_from_fasta:
                result, similar_virus_id = filter_seq_with_different_probability_of_nucleotides_from_fasta(sam_record,
                                                                                                           fasta_files,
                                                                                                           deviation_for_filtering_probability_from_fasta)
                if not result:
                    count_of_filtered_out += 1
                    continue

            # Filtering based on repeating nucleotides
            if filter_repeating and not calc_longest_repeating_subsequence(sam_record.SEQ, max_percent_limit,
                                                                           max_repeating_size):
                count_of_filtered_out += 1
                continue

            candidates.append(sam_record)

        return candidates


def decode_cigar_string(cigar):  # type: (str) -> list[str]
    """
    :param cigar: Cigar string which will be converted to string list
    :return: list of possible matches [M, M, M, M, ...], of max size of len sequence
    """
    cigar_list = []
    type_of_match = ""
    count = ""
    for ch in cigar:
        if ch.isalpha():
            type_of_match = ch

            if count != "" and type_of_match != "":
                for j in range(0, int(count)):
                    cigar_list.append(type_of_match)

                type_of_match = ""
                count = ""

        elif ch.isnumeric():
            count = count + ch

    return cigar_list



def calc_coverage_array_for_gene(count_of_genes_in_which_is_record, count_of_genes_in_which_is_not_record, gene, sam_record):  # type: (int, int, Gene, SamRecord) -> (int, int)
    """
    :param count_of_genes_in_which_is_record: variable to store the count of genes in
    :param count_of_genes_in_which_is_not_record: variable to store the count of genes out
    :param gene: Checking each record with each gene
    :param sam_record: SAM record
    :return: Count of genes in and out, and initialized coverage array [Covered, Covered, Not Covered, ...] from cigar.
    """


    # gene.record.start - len(sam_record.SEQ) -> Because mapped seq can overlap
    if gene.record.start - len(sam_record.SEQ) <= sam_record.POS < gene.record.end:
        cigar_list = decode_cigar_string(sam_record.CIGAR)
        if any(ch == "M" for ch in cigar_list):
            count_of_genes_in_which_is_record += 1
            #if not sam_record in gene.hits:
            gene.hits.append(sam_record)
            cigar_counter = 0
            for gene_pos in range(sam_record.POS, sam_record.POS + len(sam_record.SEQ)):
                if gene_pos in gene.coverage_array.keys() and cigar_list[cigar_counter] == "M":
                    gene.coverage_array[gene_pos] = True
                cigar_counter += 1
        else:
            count_of_genes_in_which_is_not_record += 1
    else:
        count_of_genes_in_which_is_not_record += 1


    return count_of_genes_in_which_is_record, count_of_genes_in_which_is_not_record


def calc_long_ends_starts(sam_record, list_of_found, n=10):  # type: (SamRecord, list[SamRecord], int) -> bool
    """
    :param sam_record: Sam record to calculate
    :param list_of_found: List of found records. New founds are appended to this list.
    :param n: Length of A^n... on start or end to look at.
    :return: True if SAM record contains A^n... on start or end and adding sam_record to list_of_found.
    """
    A = "".join(["A" for i in range(0, n)])
    C = "".join(["C" for i in range(0, n)])
    G = "".join(["G" for i in range(0, n)])
    T = "".join(["T" for i in range(0, n)])
    possibilities = (A, C, G, T)
    if sam_record.SEQ.startswith(possibilities) or sam_record.SEQ.endswith(possibilities):
        list_of_found.append(sam_record)
        return True
    else:
        return False


def filter_seq(sam_record, complex_filter_for_blast=False, deviation=0.05):  # type:(SamRecord, bool, float) -> bool
    result = False

    if complex_filter_for_blast:
        result = filter_seq_with_n(sam_record) and filter_seq_with_different_probability_of_nucleotides(sam_record,
                                                                                                        deviation)
    else:
        result = filter_seq_with_n(sam_record)

    return result


def filter_seq_with_different_probability_of_nucleotides(sam_record, deviation=0.05):  # type:(SamRecord, float) -> bool
    """
    :param sam_record: SAM record
    :param deviation: deviation to take account
    :return: True if sequence has probability of each nucleotide in the range 0.25 - deviation <= probability <= 0.25 + deviation
    """
    nucleotides_probability = {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0}
    len_of_seq = len(sam_record.SEQ)
    for nucleotide in nucleotides_probability.keys():
        count = sam_record.SEQ.count(nucleotide)
        nucleotides_probability[nucleotide] = count / len_of_seq

    if all([0.25 - deviation <= x <= 0.25 + deviation for x in nucleotides_probability.values()]):
        return True
    return False


def filter_seq_with_different_probability_of_nucleotides_from_fasta(sam_record, fasta_files, deviation=0.05):  # type:(SamRecord, list[FastaFile], float) -> (bool, list[str])
    if fasta_files is None:
        print("[SAM Analyzer]: ERROR: FASTA files are missing.")
        return False

    # probability of nucleotides in sam_record
    nucleotides_probability_sam_record = {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0}
    len_of_seq = len(sam_record.SEQ)
    for nucleotide in nucleotides_probability_sam_record.keys():
        count = sam_record.SEQ.count(nucleotide)
        nucleotides_probability_sam_record[nucleotide] = count / len_of_seq

    # probability of each fasta file
    similar_viruses_id = []
    for fasta_file in fasta_files:
        nucleotides_probability_fasta = fasta_file.nucleotides_probability
        good = True
        for nucleotide in nucleotides_probability_fasta:
            if nucleotides_probability_fasta[nucleotide] - deviation <= nucleotides_probability_sam_record[nucleotide] <= nucleotides_probability_fasta[nucleotide]:
                good = True
            else:
                good = False
        if good:
            similar_viruses_id.append(fasta_file.virus_id)

    if len(similar_viruses_id) > 0:
        return True, similar_viruses_id
    else:
        return False, similar_viruses_id


def filter_seq_with_n(sam_record):  # type:(SamRecord) -> bool
    """
    :param sam_record: SAM Record
    :return: True if seq doesn't contain "N", else False
    """
    if 'N' in sam_record.SEQ:
        return False

    return True


def calc_longest_repeating_subsequence(sequence, max_percent_limit, max_size):  # type:(str, float, int) -> bool
    """
    :param sequence: Sequence of A,C,G,T
    :param max_percent_limit: <0.0, 1.0>, Max percent limit. If percent of covered sequence by repeating seq is greater than this limit, then seq is filtered out.
    :param max_size: Size of repeating substring <2, infinity)
    :return: True if sequence is filtered out.
    """
    res, length, count = LRS(sequence, max_size)
    total_repeated_nucleotides_percent = (length * count) / len(sequence)

    if total_repeated_nucleotides_percent <= max_percent_limit:
        return True
    else:
        return False


def LRS(sequence, max_size):  # type:(str, int) -> (str, int, int)
    """
    :param sequence: Sequence of A,C,G,T
    :param max_size: Size of repeating substring <2, infinity)
    :return: Repeating substring, length of repeating substring, count of occurences of repeating string in original sequence.
    """
    str = sequence
    n = len(str)
    LCSRe = [[0 for x in range(n + 1)] for y in range(n + 1)]

    res = ""  # To store result
    res_length = 0  # To store length of result

    # building table in bottom-up manner
    index = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):

            # (j-i) > LCSRe[i-1][j-1] to remove
            # overlapping
            if (str[i - 1] == str[j - 1] and LCSRe[i - 1][j - 1] < (j - i)):

                LCSRe[i][j] = LCSRe[i - 1][j - 1] + 1
                # updating maximum length of the
                # substring and updating the finishing
                # index of the suffix
                # LCSRe[i][j] <= max_size#:
                if (LCSRe[i][j] > res_length):
                    res_length = LCSRe[i][j]
                    index = max(i, index)


            else:
                LCSRe[i][j] = 0

    # If we have non-empty result, then insert
    # all characters from first character to
    # last character of string
    if (res_length > 0):
        for i in range(index - res_length + 1, index + 1):
            res = res + str[i - 1]

    results = 0
    sub_len = len(res)
    for i in range(len(str)):
        if str[i:i + sub_len] == res:
            results += 1

    return res, len(res), str.count(res)


def create_html_output(data, filename):  # type:(DataForHTMLOutput, str) -> None
    #cmd = "pip install dominate"
    #subprocess.call(cmd, shell=False)

    shutil.copyfile(os.path.join("Core", "Classes", "functions.js"), os.path.join(data.default_output, "js", "functions.js"))
    doc = dominate.document(title='SAM file analyzer output')

    with doc.head:
        # link(rel='stylesheet', href='style.css')
        meta(charset="utf-8")
        meta(name="viewport", content="width=device-width, initial-scale=1")
        link(rel="stylesheet", href="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css")
        script(src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js")

        #
        script(src = "https://code.jquery.com/jquery-3.6.0.min.js", integrity = "sha256-/xUj+3OJU5yExlq6GSYGSHk7tPXikynS7ogEvDej/m4=", crossorigin = "anonymous")
        script(src = "https://cdnjs.cloudflare.com/ajax/libs/jqueryui/1.12.1/jquery-ui.min.js", integrity = "sha512-uto9mlQzrs59VwILcLiRYeLKPPbS/bT71da/OEBYEwcdNUk8jYIy+D176RYoop1Da+f9mvkYrmj5MCLZWEtQuA==",
               crossorigin = "anonymous", referrerpolicy = "no-referrer")
        #
        script(src="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js")

        #script(src="https://cdn.jsdelivr.net/npm/chart.js/dist/chart.min.js")
        #script(src="https://cdn.jsdelivr.net/npm/chartjs-plugin-datalabels/dist/chartjs-plugin-datalabels.min.js")
        #script(src="https://cdn.jsdelivr.net/npm/chart.js@3.0.0/dist/chart.min.js")
        #script(src="https://cdn.jsdelivr.net/npm/chartjs-plugin-datalabels@2.0.0")
        # script(src="https://cdn.jsdelivr.net/npm/bootstrap@5.0.2/dist/js/bootstrap.bundle.min.js", integrity="sha384-MrcW6ZMFYlzcLA8Nl+NtUVF0sA7MsXsP1UyJoMp4YLEuNSfAP+JcXn/tWtIaxVXM", crossorigin="anonymous")

        # link(href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css", rel="stylesheet", integrity="sha384-1BmE4kWBq78iYhFldvKuhfTAU6auU8tT94WrHftjDbrCEXSU1oBoqyl2QvZ6jIW3", crossorigin="anonymous" )
        # script(type='text/javascript', src='script.js')

    with doc.body:

        with ul(_class="nav nav-tabs"):
            # TAB WIDGET LINKS
            menu_link = li(_class="active nav-link").add(a("Summary", href="#Summary", _class="nav-link"))
            menu_link["data-toggle"] = "tab"

            menu_link = li(_class="nav-link").add(a("Virus coverage", href="#VirusCoverage", _class="nav-link"))
            menu_link["data-toggle"] = "tab"

            menu_link = li(_class="nav-link").add(a("Gene coverage", href="#GeneCoverage", _class="nav-link", id="menu-GeneCoverage"))
            menu_link["data-toggle"] = "tab"

            if data.find_sam_records_long_ends_starts:
                menu_link = li(_class="nav-link").add(
                    a("SAM records with long repeating end/start", href="#LongEndsStarts", _class="nav-link"))
                menu_link["data-toggle"] = "tab"
            if data.blast_is_set:
                menu_link = li(_class="nav-link").add(
                    a("Repeating results from blast api.", href="#BlastApi", _class="nav-link"))
                menu_link["data-toggle"] = "tab"

        # TAB WIDGET CONTENT
        with div(_class="tab-content"):
            # SUMMARY TAB
            with div(id="Summary", _class="tab-pane fade in active"):
                with div(_class="container"):
                    is_more = False
                    with ul(_class="list-group", style="margin:auto;width:60%;max-width:500px;min-width:400px;"):
                        li("Summary of SAM file records", _class="list-group-item active")
                        with li(f"Total count of records: ", _class="list-group-item"):
                            strong(data.total_count_of_sam_records)
                        with li(f"Count of records without gtf file: ", _class="list-group-item"):
                            strong(data.count_of_sam_records_with_no_gtf)
                        with li(f"Count of records containing 'N': ", _class="list-group-item"):
                            strong(data.count_of_sam_records_filtered_out_with_n)
                        with li(f"Count of records mapped on some gtf file: ", _class="list-group-item"):
                            strong(data.count_of_mapped_in_sam_records)
                        with li(f"Count of records not mapped on some gtf file: ", _class="list-group-item"):
                            strong(data.count_of_not_mapped_records)
                        if data.blast_is_set:
                            with li(f"Count of candidates from not mapped results for Blast Api: ",
                                    _class="list-group-item"):
                                strong(data.count_of_candidates_for_blast)
                                span(f", only {len(data.best_candidates)} candidates were sent to api")
                            with li(f"Count of filtered out candidates for Blast Api: ", _class="list-group-item"):
                                strong(data.count_of_filtered_out_candidates)
                        # li()
                    with table(_class="table"):
                        caption("Summary of found viruses in SAM file")
                        with thead():
                            header = tr()
                            header.add(th("#", _class="col"))
                            header.add(th("Virus ID", _class="col"))
                            header.add(th("Virus name", _class="col"))
                            header.add(th("*Count of all occurrences", _class="col"))
                            header.add(th("Count of ambiguous occurrences", _class="col"))
                            header.add(th("Count of mapped on any gtf file occurrences", _class="col"))
                            header.add(th("Count of not mapped on any gtf file occurrences", _class="col"))

                        with tbody():
                            for no, item in enumerate(data.viruses_with_counts):
                                row = tr()
                                if no > 9:
                                    row["class"] = "summary_virus_more"
                                    row["style"] = "display:none;"
                                    is_more = True
                                row.add(td(no + 1))
                                row.add(td(item[0].virus_id))
                                row.add(td(item[0].virus_name))
                                row.add(td(item[1]))
                                row.add(td(item[2]))
                                row.add(td(item[3]))
                                row.add(td(item[4]))

                    small(
                        "*(includes ambiguous occurrences, mapped and not mapped records but excludes records containing 'N')")
                    if is_more:
                        button("Show all", type="button", _class="btn btn-primary", onclick="showMore(\'summary_virus_more\')",
                               id="summary_virus_more", style="margin:3% 50% 3% 50%;")

            # VIRUS COVERAGE TABLE + GRAPHS
            with div(id='VirusCoverage', _class="tab-pane fade"):
                with div(_class="container"):
                    canvas(id="virusChart", style="width:100%;min-width:700px;max-width:800px;margin:auto;")
                    p("This graph shows only 10 viruses with heighest percent coverage value, all results are in table.",
                      style="margin: 10px;")
                    is_more = False
                    with table(_class="table"):
                        caption("Viruses and their coverage in percents of their all genes")
                        with thead():
                            header = tr()
                            header.add(th("#", _class="col"))
                            header.add(th("Virus ID", _class="col"))
                            header.add(th("Virus name", _class="col"))
                            header.add(th("Percentage of covered", _class="col"))
                            header.add(th("Hits", _class="col"))
                            header.add(th("Count of covered nucleotides", _class="col"))
                            header.add(th("Total count of nucleotides (covered+not covered)", _class="col"))
                        with tbody():
                            for i, c_v in enumerate(data.virus_coverage):
                                row = tr()
                                if i > 9:
                                    row["class"] = "virus_more"
                                    row["style"] = "display:none;"
                                    is_more = True

                                row.add(th(i + 1, scope="row"))
                                row.add(td(c_v[0]))
                                row.add(td(a(data.map_virus_id_name[c_v[0]],
                                           onclick=f"hashchanged(\"#GeneCoverage_{c_v[0]}\")", _class="tab-link")))
                                #menu_link = li(_class="nav-link").add(
                                #    a("Virus coverage", href="#VirusCoverage", _class="nav-link"))
                                #menu_link["data-toggle"] = "tab"
                                row.add(td(round(c_v[1], 2)))
                                row.add(td(len(c_v[4])))
                                row.add(td(c_v[3]))
                                row.add(td(c_v[2]))
                    if is_more:
                        button("Show all", type="button", _class="btn btn-primary", onclick="showMore(\'virus_more\')",
                               id="virus_more", style="margin:3% 50% 3% 50%;")

                    with div(_class="container"):
                        for no in range(0, len(data.virus_coverage) + 2, 2):
                            with div(_class="row"):
                                if no + 0 < len(data.virus_coverage):
                                    with div(_class="col-sm-6"):
                                        canvas(id=f"virus_{data.virus_coverage[no + 0][0]}", )
                                else:
                                    break
                                if no + 1 < len(data.virus_coverage):
                                    with div(_class="col-sm-6"):
                                        canvas(id=f"virus_{data.virus_coverage[no + 1][0]}", )
                                else:
                                    break

            # GENES COVERAGE TABLE
            with div(id='GeneCoverage', _class="tab-pane fade"):
                with div(_class="container"):
                    is_more = False
                    with table(_class="table"):
                        caption("Each gene coverage in percents")
                        with thead():
                            header = tr()
                            header.add(th("#", _class="col"))
                            header.add(th("Gene ID", _class="col"))
                            header.add(th("Protein ID", _class="col"))
                            header.add(th("Virus ID", _class="col"))
                            header.add(th("Virus name", _class="col"))
                            header.add(th("Percentage of covered", _class="col"))
                            header.add(th("Hits", _class="col"))
                            header.add(th("Count of covered nucleotides", _class="col"))
                            header.add(th("Total count of nucleotides (covered+not covered)", _class="col"))
                        with tbody():
                            c = 0
                            for virus_id, genes in data.genes_coverage.items():
                                for gene in genes:
                                    row = tr()
                                    row["class"] = ""
                                    if c > 9:
                                        row["class"] = "gene_more"
                                        row["style"] = "display:none;"
                                        is_more = True
                                    c += 1
                                    row.add(td(c, scope="row"))
                                    row.add(td(gene.record.attributes[
                                                   'gene_id'] if 'gene_id' in gene.record.attributes.keys() else ''))
                                    row.add(td(gene.record.attributes[
                                                   'protein_id'] if 'protein_id' in gene.record.attributes.keys() else ''))
                                    row.add(td(virus_id))
                                    row.add(td(data.map_virus_id_name[virus_id]))
                                    row.add(td(round(gene.covered_percents, 2)))
                                    row.add(td(len(gene.hits)))
                                    row.add(td(gene.count_of_covered_nucleotides))
                                    row.add(td(gene.length_of_gene))
                                    row["class"] += " " + virus_id
                                    row["class"] += " gene_record"

                    if is_more:
                        button("Show all", type="button", _class="btn btn-primary", onclick="showMore(\'gene_more\')",
                               id="gene_more", style="margin:3% 50% 3% 50%;")

            # SAM RECORD WITH LONG ENDS/START
            if data.find_sam_records_long_ends_starts:
                with div(id="LongEndsStarts", _class="tab-pane fade"):
                    with div(_class="container-fluid"):
                        is_more = False
                        with table(_class="table"):
                            caption("SAM records with Long repeating ends or starts")
                            with thead():
                                header = tr()
                                header.add(th("#", _class="col"))
                                header.add(th("QNAME", _class="col"))
                                header.add(th("RNAME", _class="col"))
                                header.add(th("CIGAR", _class="col"))
                                header.add(th("SEQ", _class="col"))
                            with tbody():
                                c = 0
                                for sam_record in data.sam_records_long_ends_starts:
                                    row = tr()
                                    if c > 9:
                                        row["class"] = "long_rep_more"
                                        row["style"] = "display:none;"
                                        is_more = True

                                    c += 1
                                    row.add(td(c, scope="row"))
                                    row.add(td(sam_record.QNAME))
                                    row.add(td(sam_record.virus.virus_id))
                                    row.add(td(sam_record.CIGAR))
                                    row.add(td(sam_record.SEQ))
                        if is_more:
                            button("Show all", type="button", _class="btn btn-primary",
                                   onclick="showMore(\'long_rep_more\')",
                                   id="long_rep_more", style="margin:3% 50% 3% 50%;")

            # INTERESTING RESULTS FROM BLAST API
            if data.blast_is_set:
                with div(id="BlastApi", _class="tab-pane fade"):
                    with div(_class="container"):
                        is_more = False
                        with ul(_class="list-group"):
                            li("Repeating results from Blast API", _class="list-group-item active")
                            c = 0
                            for res, count in data.interesting_results_from_api.items():
                                row = li(f"Found {count}x {res}", _class="list-group-item")
                                if c > 9:
                                    row["class"] = "list-group-item blast_res"
                                    row["style"] = "display:none;"
                                    is_more = True

                                c += 1
                        if is_more:
                            button("Show all", type="button", _class="btn btn-primary",
                                   onclick="showMore(\'blast_res\')",
                                   id="blast_res", style="margin:3% 50% 3% 50%;")


        # script(src="https://cdn.jsdelivr.net/npm/bootstrap@5.0.2/dist/js/bootstrap.bundle.min.js", integrity="sha384-MrcW6ZMFYlzcLA8Nl+NtUVF0sA7MsXsP1UyJoMp4YLEuNSfAP+JcXn/tWtIaxVXM", crossorigin="anonymous")
        script(type="text/javascript", src=os.path.join("json", "virus_coverage_output.json"))
        script(type="text/javascript", src=os.path.join("json", "genes_coverage_output.json"))

        #script(type="text/javascript", src=data.default_json_v)
        #script(type="text/javascript", src=data.defaul_json_g)
        #script(src="https://cdn.jsdelivr.net/npm/chart.js/dist/chart.min.js")
        #script(src="https://cdn.jsdelivr.net/npm/chartjs-plugin-datalabels/dist/chartjs-plugin-datalabels.min.js")
        script(src="https://cdn.jsdelivr.net/npm/chart.js@3.0.0/dist/chart.min.js")
        #script(src="https://cdn.jsdelivr.net/npm/chartjs-plugin-datalabels@2.0.0")
        #script(src=os.path.join("js", "chart.min.js"))
        #script(src=os.path.join("js", "chartjs-pligin-datalabels.js"))
        script(type="text/javascript", src=os.path.join("js", "functions.js"))

    now = datetime.datetime.now()
    file = f'{os.path.basename(filename).replace(".all.sam","")}_{now.month}-{now.day}-{now.hour}-{now.minute}.html'
    with open(os.path.join(data.default_output, file), 'w') as f:
        f.write(doc.render())
    print("[SAM Analyzer]: -----------------------------------------------")
    print(f"[SAM Analyzer]: Output was saved in file /{data.default_output}/{file}")
    print("[SAM Analyzer]: -----------------------------------------------")
