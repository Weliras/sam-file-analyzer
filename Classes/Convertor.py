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
                        gene.records.append(deepcopy(gtf_line))
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
    def load_file(url:str, map:dict[str, str], genes:list[Gene]) -> list[SamRecord]:
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
        :return: list[Virus, count_of_all, count_of_amb], grouped by attr
        """
        viruses_with_count = []
        for sam_record in sam_records:

            #if sam_record.virus.virus_name == "Human_papillomavirus_71":
            #    print()
            #    print()

            count_of_genes_in_which_is_record = 0
            count_of_genes_in_which_is_not_record = 0
            for gene in sam_record.virus.genes:
                for gene_record in gene.records:
                    if gene_record.feature == feature:
                        if gene_record.start - len(sam_record.SEQ) <= sam_record.POS < gene_record.end:
                            count_of_genes_in_which_is_record += 1
                        else:
                            count_of_genes_in_which_is_not_record += 1

            if not any(s for s in viruses_with_count if getattr(s[0], attr) == getattr(sam_record.virus, attr)):
                tmp = [sam_record.virus, 1, 0, 0, 0]
                if count_of_genes_in_which_is_record > 0:
                    tmp[3] = 1
                elif count_of_genes_in_which_is_not_record > 0:
                    tmp[4] = 1
                viruses_with_count.append(tmp)
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
                for gene in amb_virus.genes:
                    for gene_record in gene.records:
                        if gene_record.feature == feature:
                            if gene_record.start - len(sam_record.SEQ) <= sam_record.POS < gene_record.end:
                                count_of_amb_genes_in_which_is_record += 1
                            else:
                                count_of_amb_genes_in_which_is_not_record += 1
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
