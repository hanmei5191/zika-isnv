from Bio import SeqIO
import sys

fastq_set = set()

for rec in SeqIO.parse(open(sys.argv[1]), "fastq"):
    if rec.id in fastq_set:
        continue
    fastq_set.add(rec.id)
    SeqIO.write([rec], sys.stdout, "fastq")

