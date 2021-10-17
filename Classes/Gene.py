import csv
from io import StringIO


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

    def load_from_line(self, line_split):
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
    def __init__(self, virus_id, records=[]):
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
        self.records = records
        self.virus_id = virus_id

    def __del__(self):
        """self.gene.clear()
        self.transcript.clear()
        self.CDS.clear()
        self.stop_codon.clear()
        self.start_codon.clear()"""
        self.virus_id = ""
        self.records.clear()

