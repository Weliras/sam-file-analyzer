from Classes.Virus import Virus


class SamRecord:
    def __init__(self, virus=None, ambiguous_viruses=[], qname=None, flag=None, rname=None, pos=None, mapq=None, cigar=None,
                 rnext=None, pnext=None, tlen=None, seq=None, qual=None):
        self.QNAME = qname
        self.FLAG = flag
        self.RNAME = rname
        self.POS = int(pos) if str.isnumeric(pos) else pos
        self.MAPQ = mapq
        self.CIGAR = cigar
        self.RNEXT = rnext
        self.PNEXT = pnext
        self.TLEN = tlen
        self.SEQ = seq
        self.QUAL = qual

        self.virus = virus          # Virus id, Virus name, ...
        self.ambiguous_viruses = ambiguous_viruses  # List of Amb. Viruses

