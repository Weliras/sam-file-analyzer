from __future__ import annotations
import csv
import sys
import traceback
from io import StringIO
from itertools import groupby


class GTF_File_Line:
    def __init__(self, seqname=None, source=None, feature=None, start=None, end=None, score=None, strand=None, frame=None,
                 attributes=None):

        if attributes is None:
            attributes = dict()
        self.seqname = seqname          # Sequence name
        self.source = source            # Program
        self.feature = feature          # CDS, start_codon, stop_codon
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand            # forward / reverse
        self.frame = frame
        self.attributes = attributes        # ; delimited attributes

    def load_from_line(self, line_split: [str]) -> None:
        self.seqname = line_split[0] if line_split[0] != "." else None
        self.source = line_split[1] if line_split[1] != "." else None
        self.feature = line_split[2] if line_split[2] != "." else None
        self.start = int(line_split[3]) if line_split[3] != "." else None
        self.end = int(line_split[4]) if line_split[4] != "." else None
        self.score = float(line_split[5]) if line_split[5] != "." else None
        self.strand = line_split[6] if line_split[6] != "." else None
        self.frame = line_split[7] if line_split[7] != "." else None

        """for line in csv.reader(StringIO(line_split[8]), quotechar='"', delimiter=';',
                     quoting=csv.QUOTE_ALL, skipinitialspace=True):
            for attr in line:
                if 'note' in attr:
                    pass
                    #print()
                    #print()

                for line2 in csv.reader(StringIO(attr.lstrip()), delimiter=" "):
                    if len(line2) < 2:
                        continue
                    self.attributes[line2[0]] = line2[1]"""
        for attr in line_split[8].split(";"):
                if len(attr.lstrip().split(" ", 1)) != 2:
                    continue
                key = attr.lstrip().split(" ")[0]
                value = "".join(attr.lstrip().split(" ")[1:]).replace("\"", "")
                self.attributes[key] = value

class Gene:
    #def __init__(self, virus_id, GENE=None , CDS=[], start_codon=[], stop_codon=[], gene=[], transcript=[]):
    def __init__(self, virus_id: str, records=[]):
        """if GENE != None:
            self.gene = GENE.gene.copy()
            self.CDS = GENE.CDS.copy()
            self.transcript = GENE.transcript.copy()
            self.start_codon = GENE.start_codon.copy()
            self.stop_codon = GENE.stop_codon.copy()
            self.virus_id = GENE.virus_id
            return

        self.gene = gene
        self.CDS = CDS
        self.transcript = transcript
        self.start_codon = start_codon
        self.stop_codon = stop_codon

        """
        if len(records) == 1:
            self._record = records[0]
        else:
            self._record = None
        self.records = records
        self.virus_id = virus_id

        self.coverage_array = dict()
        self.length_of_gene = 0

        self.covered_percents = 0

    @property
    def record(self):
        return self._record

    @record.setter
    def record(self, value: GTF_File_Line):
        self._record = value
        self.coverage_array.clear()
        self.length_of_gene = self._record.end - self._record.start
        for i in range(self._record.start, self._record.end):                # Including upper and not Including lower bound
            self.coverage_array[i] = False


    def __del__(self):
        """self.gene.clear()
        self.transcript.clear()
        self.CDS.clear()
        self.stop_codon.clear()
        self.start_codon.clear()"""
        self.virus_id = ""
        self.records.clear()

    @staticmethod
    def write_to_file_genes_with_percents(genes: [Gene], only_non_empty: bool = False,
                                          filename: str = "gene_coverage_output.txt") -> dict[str, list[Gene]]:
        try:
            result = []
            with open(filename, "w") as file:
                for gene in genes:
                    #size_of_gene = len(gene.coverage_array.values())  # Not including end
                    size_of_gene = gene.length_of_gene
                    #if size_of_gene <= 0:
                    #    size_of_gene = 1
                    count_of_covered = sum([1 for c in gene.coverage_array.values() if c == True])
                    if only_non_empty and count_of_covered > 0:
                        file.write(
                            f"Gene Id: {gene.record.attributes['gene_id'] if 'gene_id' in gene.record.attributes.keys() else ''}"
                            f"\t Protein Id: {gene.record.attributes['protein_id'] if 'protein_id' in gene.record.attributes.keys() else ''}"
                            f"\t Virus Id: {gene.virus_id}"
                            f"\t Percentage of covered: {count_of_covered / size_of_gene * 100} %"
                            f"\t Percentage of not covered: {(size_of_gene - count_of_covered) / size_of_gene * 100} %"
                            f"\t From: {size_of_gene}"
                            f"\n")
                        result.append(gene)
                        gene.covered_percents = count_of_covered / size_of_gene * 100
                    elif not only_non_empty:
                        file.write(
                            f"Gene Id: {gene.record.attributes['gene_id'] if 'gene_id' in gene.record.attributes.keys() else ''}"
                            f"\t Protein Id: {gene.record.attributes['protein_id'] if 'protein_id' in gene.record.attributes.keys() else ''}"
                            f"\t Virus Id: {gene.virus_id}"
                            f"\t Percentage of covered: {count_of_covered / size_of_gene * 100} %"
                            f"\t Percentage of not covered: {(size_of_gene - count_of_covered) / size_of_gene * 100} %"
                            f"\t From: {size_of_gene}"
                            f"\n")
                        gene.covered_percents = count_of_covered / size_of_gene * 100
                        result.append(gene)

            # Sort by percents and group by virus id
            result.sort(key=lambda gene: (gene.covered_percents,), reverse=True)

            result = [(item, [j for j in group]) for (item, group) in groupby(result, lambda gene: gene.virus_id)]
            result2 = {}
            for v_id, genes in result:
                if v_id in result2.keys():
                    result2[v_id] += genes
                else:
                    result2[v_id] = genes
            return result2

        except Exception as e:
            print(e)
            traceback.print_exc(file=sys.stdout)
            return {}

    @staticmethod
    def write_to_file_virus_with_percents(genes: [Gene], only_non_empty: bool = False,
                                          filename: str = "virus_coverage_output.txt") -> list[list[int, float, int]]:
        try:
            # [Virus_id, %, count]
            result = []
            for gene in genes:

                #if gene.virus_id == "NC_001355.1":
                #    print(end="")
                #    print(end="")

                # count of covered pos of this gene
                count_of_covered = sum([1 for c in gene.coverage_array.values() if c == True])
                if only_non_empty and count_of_covered > 0:
                    if not any(s for s in result if s[0] == gene.virus_id):
                        size_of_gene = gene.length_of_gene
                        percentage_of_covered = (count_of_covered / size_of_gene) * 100

                        tmp = [gene.virus_id, percentage_of_covered, size_of_gene]
                        result.append(tmp)
                    else:
                        ind = [result.index(item) for item in result if item[0] == gene.virus_id]

                        # size of actual gene
                        size_of_gene = gene.length_of_gene

                        # count of already computed mapped genes
                        count_of_already_covered = (result[ind[0]][1] / 100) * result[ind[0]][2]

                        # Percentage of all covered
                        percentage_of_covered = (count_of_covered + count_of_already_covered) / (result[ind[0]][2] + size_of_gene) * 100

                        result[ind[0]][2] += size_of_gene
                        result[ind[0]][1] = percentage_of_covered

            result.sort(key=lambda virus: virus[1], reverse=True)

            with open(filename, "w") as file:
                for res in result:
                    file.write(f"Virus Id: {res[0]}\t"
                               f"Percentage of covered: {res[1]} %\t"
                               f"From: {res[2]}\n")

            # [gene.virus_id, percentage_of_covered, size_of_gene]
            return result
        except Exception as e:
            print(e)
            traceback.print_exc(file=sys.stdout)
            return []
