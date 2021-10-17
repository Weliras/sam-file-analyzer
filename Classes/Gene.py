

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
        for attr in line_split[8].split(";"):
            if len(attr.lstrip().split(" ")) != 2:
                break
            key = attr.lstrip().split(" ")[0]
            value = attr.lstrip().split(" ")[1]
            self.attributes[key] = value


class Gene:
    def __init__(self, virus_id, CDS=None, start_codon=None, stop_codon=None):
        self.CDS = CDS
        self.start_codon = start_codon
        self.stop_codon = stop_codon

        self.virus_id = virus_id


