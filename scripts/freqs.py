import pysam
import sys
import csv 

def read_variants(fn):
    #Name    Minimum Maximum Length  Change  Coverage    Polymorphism Type   Variant Frequency   replica modality    freq
    snps = {}

    for row in csv.DictReader(open(fn), dialect='excel-tab'):
        if row['Polymorphism Type'].startswith('SNP'):
            a, b = row['Change'].split(" -> ")

            snps[int(row['Minimum'])] = (a, b)

    return snps

snps = read_variants(sys.argv[1])

print ("Pos\tQual\tFreq\tBase")

samfile = pysam.AlignmentFile(sys.argv[2], "rb")
for pileupcolumn in samfile.pileup(sys.argv[3]):
    for q in [0, 5, 10,11,12,13,14,15,16,17,18,19,20]:
        freqs = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        if pileupcolumn.pos+1 in snps:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    # query position is None if is_del or is_refskip is set.

                    if pileupread.alignment.query_qualities[pileupread.query_position] >= q:
                         freqs[pileupread.alignment.query_sequence[pileupread.query_position]] += 1
            if max(freqs.values()):
                nonindel_coverage = sum(freqs.values())

                snp = snps[pileupcolumn.pos+1]
                if snp[1] in ('A', 'T', 'G', 'C') :
                    # ignore multi-allelic sites
                #print ("%s\t%s\t%s\t%s" % (pileupcolumn.pos+1, q, float(max(freqs.values())) / (nonindel_coverage)), )
                    print ("%s\t%s\t%s\t%s" % (pileupcolumn.pos+1, q, float(freqs[snp[1]] / nonindel_coverage), snp[1]))


