import os
import urllib.request
import sys, traceback
from copy import deepcopy

from Classes.Gene import Gene, GTF_File_Line
from Classes.Virus import Virus
from Classes.SamRecord import SamRecord


class Convertor:

    @staticmethod
    def get_map_of_id_to_name(url:str) -> dict[str, str]:
        try:
            virusSeq_to_name_map = {}
            with urllib.request.urlopen(url) as file:
                for line in file:
                    decoded_line = line.decode("utf-8")
                    decoded_line = decoded_line.strip()
                    line_split = decoded_line.split('\t')
                    virusSeq_to_name_map[line_split[0]] = line_split[1]
                    # print(decoded_line)
                virusSeq_to_name_map['*'] = '*'  # For unknown
        except Exception as e:
            print(e)
            traceback.print_exc(file=sys.stdout)
        return virusSeq_to_name_map

    @staticmethod
    def load_gtf_files_only_cds_gene() -> list[Gene]:
        '''
        Function returns only cds and gene lines from gtf files
        :return: list of genes
        '''
        d = 'gtf_files'
        folders = [os.path.join(d, o) for o in os.listdir(d) if os.path.isdir(os.path.join(d, o))]
        genes = []
        for dir in folders:
            virus_id = dir.split("\\")[1].replace("-", "|")
            gtf_file = [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f))]
            if len(gtf_file) <= 0:
                print(f"Gtf file for virus_id: {virus_id} was not found.")
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

                        if gtf_line.feature != "gene" and gtf_line.feature != "CDS":
                            continue
                        gene = Gene(virus_id)
                        gene.records.append(deepcopy(gtf_line))         # deprecated
                        gene.record = deepcopy(gtf_line)
                        genes.append(deepcopy(gene))

                        del gene
            except Exception as e:
                print(e)
                traceback.print_exc(file=sys.stdout)
        return genes

    @staticmethod
    def load_gtf_files() -> list[Gene]:
        d = 'gtf_files'
        folders = [os.path.join(d, o) for o in os.listdir(d) if os.path.isdir(os.path.join(d, o))]

        genes = []
        for dir in folders:
            virus_id = dir.split("\\")[1].replace("-", "|")
            gtf_file = [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f))]
            if len(gtf_file) <= 0:
                print(f"Gtf file for virus_id: {virus_id} was not found.")
                continue
            gtf_file = gtf_file[0]
            try:
                with open(os.path.join(dir, gtf_file), "r") as file:
                    prev_protein_id = ""
                    prev_gene_id = ""
                    first_in_file = True
                    gene = Gene(virus_id)
                    for line in file:

                        if line == "###\n":                 # End of file, insert gene and delete
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
                print(e)
                print(virus_id)
                traceback.print_exc(file=sys.stdout)
        return genes

    @staticmethod
    def load_sam_files(url:str, map:dict[str, str], genes:list[Gene]) -> list[SamRecord]:
        """
        Using dynamic programming to get runtime under 30 secs
        :param genes: Loaded genes from GTF files
        :param url: Path to input SAM file
        :param map: Map of Virus id -> Virus Name
        :return: list of SamRecords loaded from a file. SamRecord->(Virus->Genes,AmbViruses->Genes)
        """
        try:
            dynamic_programming = dict()
            for virus_id in map.keys():
                dynamic_programming[virus_id] = None

            with urllib.request.urlopen(url) as file:
                sam_records = []
                for line in file:
                    decoded_line = line.decode("utf-8")
                    line_split = decoded_line.split('\t')

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
                                      ambiguous_viruses=ambiguous_viruses, pos=line_split[3], cigar=line_split[5], seq=line_split[9]))
                        continue
                    g = [g for g in genes if g.virus_id == line_split[2]]
                    sam_records.append(
                        SamRecord(virus=Virus(virus_id=line_split[2], virus_name=map[line_split[2]], genes=g),
                                  ambiguous_viruses=ambiguous_viruses,  pos=line_split[3], cigar=line_split[5], seq=line_split[9]))
                    dynamic_programming[line_split[2]] = g
        except Exception as e:
            print(e)
            traceback.print_exc(file=sys.stdout)
        return sam_records


    @staticmethod
    def get_seqs_with_count_grouped_by(sam_records:list[SamRecord], attr:str, feature="CDS") -> list[Virus, int, int, int, int]:
        """
        :param feature: "CDS" or "gene"
        :param sam_records: Loaded lines by method load_file()
        :param attr: Attribute of Virus class which is used to group viruses
        :return: list[Virus, count_of_all, count_of_amb, count_of_mapped_id, count_of_mapped_out], grouped by attr
        """
        viruses_with_count = []
        for sam_record in sam_records:

            #if sam_record.virus.virus_name == "Human_papillomavirus_71":
            #    print()
            #    print()

            # Filtering based on repeating nucleotides
            if not calc_longest_repeating_subsequence(sam_record.SEQ, 0.30):
                pass

            count_of_genes_in_which_is_record = 0
            count_of_genes_in_which_is_not_record = 0

            # Decoding CIGAR string and computing count of nucleotides which were mapped in and out.
            for gene in sam_record.virus.genes:
                if gene.record.feature == feature:
                    count_of_genes_in_which_is_record, count_of_genes_in_which_is_not_record = \
                        calc_coverage_array_for_gene(count_of_genes_in_which_is_record,
                                                              count_of_genes_in_which_is_not_record, gene=gene,
                                                              sam_record=sam_record, attr=attr)

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
                    viruses_with_count[ind[0]][3] += 1          # Add count of IN
                elif count_of_genes_in_which_is_not_record > 0:
                    viruses_with_count[ind[0]][4] += 1          # Add count of OUT
                viruses_with_count[ind[0]][1] += 1


            # counting occurences of virus in other amb viruses
            for amb_virus in sam_record.ambiguous_viruses:
                count_of_amb_genes_in_which_is_record = 0
                count_of_amb_genes_in_which_is_not_record = 0

                # Calculating coverage array of amb viruses.
                for gene in amb_virus.genes:
                    if gene.record.feature == feature:
                        count_of_amb_genes_in_which_is_record, count_of_amb_genes_in_which_is_not_record =\
                            calc_coverage_array_for_gene(count_of_amb_genes_in_which_is_record
                                                                   , count_of_amb_genes_in_which_is_not_record,
                                                                   gene=gene, attr=attr, sam_record=sam_record)

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
                    viruses_with_count[ind[0]][1] += 1      # count of all
                    viruses_with_count[ind[0]][2] += 1      # count of amb
                    if count_of_amb_genes_in_which_is_record > 0:
                        viruses_with_count[ind[0]][3] += 1      # count of amb IN
                    elif count_of_amb_genes_in_which_is_not_record > 0:
                        viruses_with_count[ind[0]][4] += 1      # count of amb OUT

        return viruses_with_count


def decode_cigar_string(cigar: str) -> [str]:
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


def calc_coverage_array_for_gene(count_of_genes_in_which_is_record:int, count_of_genes_in_which_is_not_record:int, gene:Gene, sam_record:SamRecord, attr:str):
    if gene.record.start - len(sam_record.SEQ) <= sam_record.POS < gene.record.end:
        cigar_list = decode_cigar_string(sam_record.CIGAR)
        if any(ch == "M" for ch in cigar_list):
            count_of_genes_in_which_is_record += 1
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

def calc_longest_repeating_subsequence(sequence:str, max_percent_limit:float) -> bool:
    res, length, count = LRS(sequence)
    total_repeated_nucleotides_percent = (length * count) / len(sequence)

    #print(sequence, res, length, count, f"{total_repeated_nucleotides_percent} %")
    if total_repeated_nucleotides_percent <= max_percent_limit:
        return True
    else:
        return False


def LRS(sequence:str, max_size:int) -> (str, int, int):
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
                if (LCSRe[i][j] > res_length) and LCSRe[i][j] <= 3:
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

