conda create -n zika-isnv
source activate zika-isnv
conda install porechop bwa
git clone git@github.com:artic-network/fieldbioinformatics.git
wget http://nanopore.s3.climb.ac.uk/Zika_iSNV_R9.4_albacore2.3.1.tgz
wget http://nanopore.s3.climb.ac.uk/Zika_iSNV_R9.5_albacore2.3.1.tgz
tar xvfz Zika_iSNV_R9.4_albacore2.3.1.tgz
cat zika-isnv/fast5_r94/workspace/pass/*.fastq > fast5_r94.pass.fastq
python scripts/dedup.py fast5_r94.pass.fastq > fast5_r94.dedup.pass.fastq
porechop -b fast5_r94.dedup.pass_porechop -t 4 -i fast5_r94.dedup.pass.fastq
bwa index refs/ZIKV_REF.fasta
bwa mem -x ont2d -t 4 refs/ZIKV_REF.fasta fast5_r94.dedup.pass_porechop/BC01.fastq | samtools view -bS - | samtools sort -o BC01.sorted.bam
bwa mem -x ont2d -t 4 refs/ZIKV_REF.fasta fast5_r94.dedup.pass_porechop/BC02.fastq | samtools view -bS - | samtools sort -o BC02.sorted.bam
bwa mem -x ont2d -t 4 refs/ZIKV_REF.fasta fast5_r94.dedup.pass_porechop/BC03.fastq | samtools view -bS - | samtools sort -o BC03.sorted.bam
samtools index BC01.sorted.bam
samtools index BC02.sorted.bam
samtools index BC03.sorted.bam
align_trim --normalise 1000 refs/ZikaAsian.scheme.bed <BC01.sorted.bam 2>/dev/null | samtools view -bS - | samtools sort - -o BC01.trimmed.sorted.bam
align_trim --normalise 1000 refs/ZikaAsian.scheme.bed <BC02.sorted.bam 2>/dev/null | samtools view -bS - | samtools sort - -o BC02.trimmed.sorted.bam
align_trim --normalise 1000 refs/ZikaAsian.scheme.bed <BC03.sorted.bam 2>/dev/null | samtools view -bS - | samtools sort - -o BC03.trimmed.sorted.bam
samtools index BC01.trimmed.sorted.bam
samtools index BC02.trimmed.sorted.bam
samtools index BC03.trimmed.sorted.bam
python scripts/freqs.py 2018.01.15_ZIKV-iSNV/variants.tsv BC01.trimmed.sorted.bam KU926309 > BC01.freqs.txt
python scripts/freqs.py 2018.01.15_ZIKV-iSNV/variants.tsv BC02.trimmed.sorted.bam KU926309 > BC02.freqs.txt
python scripts/freqs.py 2018.01.15_ZIKV-iSNV/variants.tsv BC03.trimmed.sorted.bam KU926309 > BC03.freqs.txt
