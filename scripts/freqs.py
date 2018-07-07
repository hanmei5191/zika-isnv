import pysam
import sys
import csv 
from Bio import SeqIO

def read_variants(fn):
    #Name    Minimum Maximum Length  Change  Coverage    Polymorphism Type   Variant Frequency   replica modality    freq
    snps = {}

    for row in csv.DictReader(open(fn), dialect='excel-tab'):
        if row['Polymorphism Type'].startswith('SNP'):
            a, b = row['Change'].split(" -> ")

            snps[int(row['Minimum'])] = (a, b)

    return snps

def go(args):
    if args.variants:
        snps = read_variants(args.variants)
    else:
        snps = None

    references = SeqIO.to_dict(SeqIO.parse(open(args.reference), "fasta"))

    print ("Pos\tQual\tFreq\tBase\tCoverage")

    samfile = pysam.AlignmentFile(args.alignment, "rb")
    for contig in references.keys():
        for pileupcolumn in samfile.pileup(contig):
            for q in [0, 5, 10,11,12,13,14,15,16,17,18,19,20]:
                freqs = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
                if not snps or (pileupcolumn.pos+1 in snps):
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            # query position is None if is_del or is_refskip is set.

                            if pileupread.alignment.query_qualities[pileupread.query_position] >= q:
                                 freqs[pileupread.alignment.query_sequence[pileupread.query_position]] += 1
                    if max(freqs.values()):
                        nonindel_coverage = sum(freqs.values())

                        if snps:
                            base = snps[pileupcolumn.pos+1][1]
                        else:
                            base = references[pileupcolumn.reference_name][pileupcolumn.reference_pos]
                            
                        if base in ('A', 'T', 'G', 'C') :
                            print ("%s\t%s\t%s\t%s\t%s" % (pileupcolumn.pos+1, q, float(freqs[base] / nonindel_coverage), base, nonindel_coverage))

import argparse

parser = argparse.ArgumentParser(description='Retrieve frequencies.')
parser.add_argument('--variants', help='Retrieve specific variants from tsv file.')
parser.add_argument('alignment', help='BAM file')
parser.add_argument('reference', help='Reference FASTA')

args = parser.parse_args()
go(args)

