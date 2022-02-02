class FastaFile:
    def __init__(self, virus_id: str, virus_name: str, sequence: str, nucleotides_probability: {str: float}):
        self.virus_id = virus_id
        self.virus_name = virus_name
        self.sequence = sequence
        self.nucleotides_probability = nucleotides_probability
