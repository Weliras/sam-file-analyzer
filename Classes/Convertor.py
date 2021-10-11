import urllib.request
import sys, traceback
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
    def load_file(url, map):
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
                                    ambiguous_viruses.append(Virus(virus_id=v[0], virus_name=map[v[0]]))
                    sam_records.append(SamRecord(virus=Virus(virus_id=line_split[2], virus_name=map[line_split[2]]),
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

