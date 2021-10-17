import os
import urllib.request
import sys, traceback

from Classes.Gene import Gene, GTF_File_Line
from Classes.Virus import Virus
from Classes.SamRecord import SamRecord


class Convertor:

    @staticmethod
    def get_map_of_id_to_name(url):
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
    def load_gtf_files():
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
                    prev_line = ""
                    gene = Gene(virus_id)
                    for line in file:
                        if line.startswith("#"):
                            continue
                        if line.split("\t")[2] == "CDS":
                            if prev_line == "CDS" or prev_line == "start_codon":
                                genes.append(Gene(virus_id=gene.virus_id, CDS=gene.CDS, start_codon=gene.start_codon, stop_codon=gene.stop_codon))

                            cds_line = GTF_File_Line()
                            line_split = line.split("\t")
                            cds_line.load_from_line(line_split)

                            gene.CDS = cds_line

                            prev_line = "CDS"
                        elif line.split("\t")[2] == "start_codon":

                            start_codon_line = GTF_File_Line()
                            line_split = line.split("\t")
                            start_codon_line.load_from_line(line_split)

                            gene.start_codon = start_codon_line

                            prev_line = "start_codon"
                        elif line.split("\t")[2] == "stop_codon":

                            stop_codon_line = GTF_File_Line()
                            line_split = line.split("\t")
                            stop_codon_line.load_from_line(line_split)

                            gene.stop_codon = stop_codon_line

                            genes.append(Gene(virus_id=gene.virus_id, CDS=gene.CDS, start_codon=gene.start_codon, stop_codon=gene.stop_codon))

                            prev_line = "stop_codon"
            except Exception as e:
                print(e)
                traceback.print_exc(file=sys.stdout)
        return genes


    @staticmethod
    def load_file(url, map, genes):
        """
        :param url: Path to input SAM file
        :param map: Map of Virus id -> Virus Name
        :return: list of SamRecords loaded from a file
        """
        try:
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
                                    # genes to virus
                                    g = [g for g in genes if g.virus_id == v[0]]
                                    ambiguous_viruses.append(Virus(virus_id=v[0], virus_name=map[v[0]], genes=g))
                    g = [g for g in genes if g.virus_id == line_split[2]]
                    sam_records.append(SamRecord(virus=Virus(virus_id=line_split[2], virus_name=map[line_split[2]], genes=g),
                                                 ambiguous_viruses=ambiguous_viruses))
        except Exception as e:
            print(e)
            traceback.print_exc(file=sys.stdout)
        return sam_records

    @staticmethod
    def get_seqs_with_count_grouped_by(sam_records, attr):
        """
        :param sam_records: Loaded lines by method load_file()
        :param attr: Attribute of Virus class which is used to group viruses
        :return: list[Virus, count_of_all, count_of_amb], grouped by attr
        """
        viruses_with_count = []
        for sam_record in sam_records:
            if not any(s for s in viruses_with_count if getattr(s[0], attr) == getattr(sam_record.virus, attr)):
                viruses_with_count.append([sam_record.virus, 1, 0])
            else:
                ind = [viruses_with_count.index(item) for item in viruses_with_count if getattr(item[0], attr) ==
                       getattr(sam_record.virus, attr)]
                viruses_with_count[ind[0]][1] += 1

            for amb_virus in sam_record.ambiguous_viruses:
                if not any(s for s in viruses_with_count if getattr(s[0], attr) == getattr(amb_virus, attr)):
                    viruses_with_count.append([amb_virus, 1, 1])
                else:
                    ind = [viruses_with_count.index(item) for item in viruses_with_count if getattr(item[0], attr) ==
                           getattr(amb_virus, attr)]
                    viruses_with_count[ind[0]][1] += 1
                    viruses_with_count[ind[0]][2] += 1

        return viruses_with_count

