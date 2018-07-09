import pysam
import sys
import csv 
from Bio import SeqIO
from copy import copy
from operator import attrgetter

class Position:
    def __init__(self, base):
        self.base = base
        self.count = 0
        self.forward = 0
        self.reverse = 0

    def update(self, is_reverse):
        self.count += 1
        if is_reverse:
            self.reverse += 1
        else:
            self.forward += 1

    def __repr__(self):
        return "%s : %s" % (self.base, self.count)

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

    print ("Pos\tQual\tFreq\tRef\tBase\tUngappedCoverage\tTotalCoverage\tVariantCov\tForwardVariantCov\tReverseVariantCov")

    samfile = pysam.AlignmentFile(args.alignment, "rb")
    for contig in references.keys():
        for pileupcolumn in samfile.pileup(contig, stepper='samtools', fastafile=fasta):
            for q in [0, ]:
                #5, 10,11,12,13,14,15,16,17,18,19,20]:
                #freqs = copy({'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0, 'R' : 0})
                freqs = {'A': Position('A'), 'T': Position('T'), 'G': Position('G'), 'C': Position('C'), '-': Position('-'), 'R': Position('R')}

                if not snps or (pileupcolumn.pos+1 in snps):
                    for pileupread in pileupcolumn.pileups:
                        if pileupread.is_del:
                            freqs['-'].update(pileupread.alignment.is_reverse)
                        elif pileupread.is_refskip:
                            freqs['R'].update(pileupread.alignment.is_reverse)
                        else:
                            # query position is None if is_del or is_refskip is set.

                            if pileupread.alignment.query_qualities[pileupread.query_position] >= q:
                                 freqs[pileupread.alignment.query_sequence[pileupread.query_position]].update(pileupread.alignment.is_reverse)

                    total_nonindel_coverage = sum((x.count for x in freqs.values()))
                    nonindel_coverage = freqs['A'].count + freqs['T'].count + freqs['G'].count + freqs['C'].count
                    if nonindel_coverage:
                        reference_name = samfile.getrname(pileupcolumn.reference_id)

                        ref = references[reference_name][pileupcolumn.reference_pos]
                        if args.snpfreqmin:
                            base = None
                            for freq in sorted(freqs.values(), reverse=True, key=attrgetter('count')):
                                if freq.base == ref:
                                    continue

                                if freq.base not in ['A', 'T', 'G', 'C']:
                                    continue

                                cov = float(freq.count) / float(nonindel_coverage)
                                if cov < args.snpfreqmin:
                                    continue
                                base = freq.base
                                break
                            if not base:
                                continue
                        elif snps:
                            base = snps[pileupcolumn.pos+1][1]
                        else:
                            base = ref

                        if base in ('A', 'T', 'G', 'C') :
                            print ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (pileupcolumn.pos+1, q, float(freqs[base].count) / float(nonindel_coverage), ref, base, nonindel_coverage, total_nonindel_coverage, freqs[base].count, freqs[base].forward, freqs[base].reverse))

import argparse

parser = argparse.ArgumentParser(description='Retrieve frequencies.')
parser.add_argument('--variants', help='Retrieve specific variants from tsv file.')
parser.add_argument('--snpfreqmin', type=float)
parser.add_argument('alignment', help='BAM file')
parser.add_argument('reference', help='Reference FASTA')

args = parser.parse_args()
go(args)

