
import requests
from Bio.Blast.NCBIWWW import qblast


class BlastApi:
    @staticmethod
    def connect():
        program = "blastn"
        database = "nt"
        seq = "AGGATATTGTATTAGACCTGCAACCTCCAGACCCTGTAGGGTTACATTGCTATGAGCAATTAGTAGACAGCGCAGA"
        sequence_data = open("blast_example.fasta").read()

        result_handle = qblast(program=program, database=database, sequence=seq, format_type="XML")

        with open('results.xml', 'w') as save_file:
            blast_results = result_handle.read()
            save_file.write(blast_results)









