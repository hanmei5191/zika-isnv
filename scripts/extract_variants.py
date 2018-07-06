import pysam
import sys
import csv
from Bio import SeqIO

class Position:
    def __init__(self, base):
        self.base = base
        self.count = 0

    def update(self):
        self.count += 1

    def __cmp__(a, b):
        return cmp(a.count, b.count)

    def __repr__(self):
        return "%s : %s" % (self.base, self.count)

QUAL_THRESHOLD = int(sys.argv[3])
COVERAGE_THRESHOLD = float(sys.argv[4])

def read_positions(fn):
	for row in csv.DictReader(open(fn), dialect='excel'):
		print row

def main():
	print """##fileformat=VCFv4.2
	##INFO=<ID=Frequency,Number=1,Type=Float,Description="The fraction of base-space reads that support the variant">
	##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample"""

	references = SeqIO.to_dict(SeqIO.parse(sys.argv[2], "fasta"))
	samfile = pysam.AlignmentFile(sys.argv[1], "rb")

	if len(sys.argv) == 6:
		positions = read_positions(sys.argv[5])
	else:
		positions = None

	for pileupcolumn in samfile.pileup():
	#    print ("\ncoverage at base %s = %s" %
	#           (pileupcolumn.pos, pileupcolumn.n))
	    freqs = {'A': Position('A'), 'T': Position('T'), 'G': Position('G'), 'C': Position('C')}


	    for pileupread in pileupcolumn.pileups:
		if not pileupread.is_del and not pileupread.is_refskip:
		    # query position is None if is_del or is_refskip is set.

		    if pileupread.alignment.query_qualities[pileupread.query_position] >= QUAL_THRESHOLD:
			freqs[pileupread.alignment.query_sequence[pileupread.query_position]].update()
	 

	    sorted_freqs = sorted(freqs.values(), reverse=True)
	    nonindel_coverage = sum([ x.count for x in sorted_freqs ])
	    first_allele = sorted_freqs[0]
	    if not nonindel_coverage:
		continue
	    first_allele_cov = float(first_allele.count) / float(nonindel_coverage)
	    if references[pileupread.alignment.reference_name][pileupcolumn.pos] != first_allele.base:
		print "%s\t%s\t%s\t%s\t%s\t%s\t%s\tFrequency=%f\t%s\t%s" % (pileupread.alignment.reference_name, pileupcolumn.pos+1, '.', references[pileupread.alignment.reference_name][pileupcolumn.pos],  first_allele.base, 100, 'PASS', first_allele_cov, 'GT', '1')
	    second_allele = sorted_freqs[1]
	    second_allele_cov = float(second_allele.count) / float(nonindel_coverage)
	    if references[pileupread.alignment.reference_name][pileupcolumn.pos] != second_allele.base and second_allele_cov >= COVERAGE_THRESHOLD:
		print "%s\t%s\t%s\t%s\t%s\t%s\t%s\tFrequency=%f\t%s\t%s" % (pileupread.alignment.reference_name, pileupcolumn.pos+1, '.', references[pileupread.alignment.reference_name][pileupcolumn.pos],  second_allele.base, 100, 'PASS', second_allele_cov, 'GT', '1')

main()

