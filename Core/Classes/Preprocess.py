import os
import sys

import subprocess
import traceback


"""
Predpokladame, ze na vstupu mame fastq(gz) soubor.
1) pokud gz, pak rozbalime
2) zalignujeme fastq vuci virove referencni skevenci bwa mem
3) profiltrujeme vysledny sam od nenamapovanych sekvenci
"""




def align(fq_path, viral_genome):
    sam_path = fq_path.replace(".fastq", ".sam")
    #print("sam_path: " + sam_path)
    #print("fq_path: " + fq_path)
    if os.path.exists(fq_path):
        cmd = "bwa mem -t 32 %s %s > %s" % (viral_genome , fq_path, sam_path)
        subprocess.call(cmd, shell=True)
    return sam_path


def filterUnmapped(sam_path):
    sam_result_path = sam_path.replace(".sam", ".all.sam")
    if os.path.exists(sam_path):
        cmd = "samtools view -F 0x4 -o %s %s" % (sam_result_path, sam_path)
        subprocess.call(cmd, shell=True)
        os.remove(sam_path)
    return sam_result_path


def createReferenceGenome(viral_genomes_folder: str, DEFAULT_REF_GENOME: str) -> str:

    all_genomes_file = DEFAULT_REF_GENOME
    if os.path.exists(viral_genomes_folder):
        folders = [os.path.join(viral_genomes_folder, o) for o in os.listdir(viral_genomes_folder) if os.path.isdir(os.path.join(viral_genomes_folder, o))]
        try:

            if os.path.exists(all_genomes_file):
                os.remove(all_genomes_file)

            with open(all_genomes_file, "a") as all_genomes:
                for dir in folders:
                    #print(dir)
                    fasta_file = [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f))]
                    if len(fasta_file) <= 0:
                        continue
                    fasta_file = fasta_file[0]
                    with open(os.path.join(dir, fasta_file), "r") as genome:
                        all_genomes.writelines(genome.readlines())
                        all_genomes.write("\n")

        except Exception as e:
            print(f"[SAM Analyzer]: {e}")
            traceback.print_exc(file=sys.stdout)

    cmd = "bwa index %s" % DEFAULT_REF_GENOME
    subprocess.call(cmd, shell=True)

    return all_genomes_file
"""
Na vstupu: 
 - adresar s fastq soubory.
 - adresar s genomy
 - viral_ids soubor
"""


def preprocess(fastq_file: str, DEFAULT_FASTA_FOLDER: str, DEFAULT_REF_GENOME: str, recreate_ref_genome: bool = True) -> str:

    print("Processing:")
    print(fastq_file)

    if os.path.exists(fastq_file) and fastq_file.endswith(".gz"):
        cmd = "gunzip " + fastq_file
        subprocess.call(cmd, shell=True)

        fastq_file = fastq_file.replace(".gz", "")

        print("[SAM Analyzer]: Fastq file gunzipped.")


    if fastq_file.endswith(".fastq"):
        fpath = fastq_file

        viral_genome_all = DEFAULT_REF_GENOME

        if recreate_ref_genome:
            viral_genome_all = createReferenceGenome(DEFAULT_FASTA_FOLDER, DEFAULT_REF_GENOME)

        #print(viral_genome_all, fpath)
        sam = align(fpath, viral_genome_all)
        mapped_sam = filterUnmapped(sam)
        return mapped_sam
