from Classes.Virus import Virus


class SamRecord:
    def __init__(self, virus, ambiguous_viruses):
        self.virus = virus          # Virus id, Virus name, ...
        self.ambiguous_viruses = ambiguous_viruses  # List of Amb. Viruses

