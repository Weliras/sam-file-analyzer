import os

from Classes.Virus import Virus
from Classes.SamRecord import SamRecord
from Classes.Convertor import Convertor

RNA_SAMPLE = "http://homel.vsb.cz/~vas218/files/viruses/viruses.filtered.sam"
VIRUS_NAME = "http://homel.vsb.cz/~vas218/files/viruses/viral_id.csv"

if __name__ == '__main__':

    # Load gtf files
    genes = Convertor.load_gtf_files()

    # Creating map of VirusSeq -> VirusName
    map = Convertor.get_map_of_id_to_name(VIRUS_NAME)

    # Load a file and mapping ids to names, list of Sequences
    # [VirusSeq, Count] -> [VirusName, Count] | VirusSeq : VirusName = N : 1
    sam_records = Convertor.load_file(RNA_SAMPLE, map, genes)

    # for sam_record in sam_records:
    #    sam_record.virus.genes = [gene for gene in genes if gene.virus_id == sam_record.virus.virus_id]

    # check for same virus names but different virus ids
    # sam_records.append(SamRecord(Virus("NC_002023.1", "Influenza_A_virus")))
    # sam_records.append(SamRecord(Virus("NC_002021.1", "Influenza_A_virus")))

    # getting list [ Virus, count_of_all, count_of_amb ] by attribute of Sequence
    virus_with_count = Convertor.get_seqs_with_count_grouped_by(sam_records, "virus_name")

    # sort it by count desc
    virus_with_count.sort(key= lambda virus: virus[1], reverse=True)

    # Printing result VirusId\tVirusName\tCountOfAll\tCountOfAmb
    print("virus_id\tvirus_name\tcount_of_all\tcount_of_amb\n")
    for item in virus_with_count:
        print(f"{item[0].virus_id}\t{item[0].virus_name}\t{item[1]}\t{item[2]}")

