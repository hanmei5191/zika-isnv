import pysam
import sys
import csv 
from Bio import SeqIO
from copy import copy

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
    fasta = pysam.Fastafile(args.reference)

    print ("Pos\tQual\tFreq\tRef\tBase\tUngappedCoverage\tTotalCoverage")

    samfile = pysam.AlignmentFile(args.alignment, "rb")
    for contig in references.keys():
        for pileupcolumn in samfile.pileup(contig, stepper='samtools', fastafile=fasta):
            for q in [0, ]:
                #5, 10,11,12,13,14,15,16,17,18,19,20]:
                freqs = copy({'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0, 'R' : 0})
                if not snps or (pileupcolumn.pos+1 in snps):
                    for pileupread in pileupcolumn.pileups:
                        if pileupread.is_del:
                            freqs['-'] += 1
                        elif pileupread.is_refskip:
                            freqs['R'] += 1
                        else:
                            # query position is None if is_del or is_refskip is set.

                            if pileupread.alignment.query_qualities[pileupread.query_position] >= q:
                                 freqs[pileupread.alignment.query_sequence[pileupread.query_position]] += 1

                    total_nonindel_coverage = sum(freqs.values())
                    nonindel_coverage = freqs['A'] + freqs['T'] + freqs['G'] + freqs['C']
                    if nonindel_coverage:
                        reference_name = samfile.getrname(pileupcolumn.reference_id)

                        ref = references[reference_name][pileupcolumn.reference_pos]
                        if snps:
                            base = snps[pileupcolumn.pos+1][1]
                        else:
                            base = ref

                        if base in ('A', 'T', 'G', 'C') :
                            print ("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (pileupcolumn.pos+1, q, float(freqs[base]) / float(nonindel_coverage), ref, base, nonindel_coverage, total_nonindel_coverage))

import argparse

parser = argparse.ArgumentParser(description='Retrieve frequencies.')
parser.add_argument('--variants', help='Retrieve specific variants from tsv file.')
parser.add_argument('alignment', help='BAM file')
parser.add_argument('reference', help='Reference FASTA')

args = parser.parse_args()
go(args)

