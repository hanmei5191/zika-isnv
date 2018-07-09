THREADS=12

source activate zika-isnv

cat zika-isnv/fast5_r94/workspace/pass/*.fastq > fast5_r94.pass.fastq
python scripts/dedup.py fast5_r94.pass.fastq > fast5_r94.dedup.pass.fastq
porechop -b fast5_r94.dedup.pass_porechop -t $THREADS -i fast5_r94.dedup.pass.fastq
bwa index refs/ZIKV_REF.fasta
bwa mem -x ont2d -t $THREADS refs/ZIKV_REF.fasta fast5_r94.dedup.pass_porechop/BC01.fastq | samtools view -bS - | samtools sort -o BC01.sorted.bam
bwa mem -x ont2d -t $THREADS refs/ZIKV_REF.fasta fast5_r94.dedup.pass_porechop/BC02.fastq | samtools view -bS - | samtools sort -o BC02.sorted.bam
bwa mem -x ont2d -t $THREADS refs/ZIKV_REF.fasta fast5_r94.dedup.pass_porechop/BC03.fastq | samtools view -bS - | samtools sort -o BC03.sorted.bam
samtools index BC01.sorted.bam
samtools index BC02.sorted.bam
samtools index BC03.sorted.bam
align_trim --normalise 1000 refs/ZikaAsian.scheme.bed <BC01.sorted.bam 2>/dev/null | samtools view -bS - | samtools sort - -o BC01.trimmed.sorted.bam
align_trim --normalise 1000 refs/ZikaAsian.scheme.bed <BC02.sorted.bam 2>/dev/null | samtools view -bS - | samtools sort - -o BC02.trimmed.sorted.bam
align_trim --normalise 1000 refs/ZikaAsian.scheme.bed <BC03.sorted.bam 2>/dev/null | samtools view -bS - | samtools sort - -o BC03.trimmed.sorted.bam
samtools index BC01.trimmed.sorted.bam
samtools index BC02.trimmed.sorted.bam
samtools index BC03.trimmed.sorted.bam
python scripts/freqs.py BC01.trimmed.sorted.bam refs/ZIKV_REF.fasta > BC01.freqs.txt
python scripts/freqs.py BC02.trimmed.sorted.bam refs/ZIKV_REF.fasta > BC02.freqs.txt
python scripts/freqs.py BC03.trimmed.sorted.bam refs/ZIKV_REF.fasta > BC03.freqs.txt
python scripts/freqs.py --variants illumina_variants.tsv BC01.trimmed.sorted.bam refs/ZIKV_REF.fasta > BC01.variants.freqs.txt
python scripts/freqs.py --variants illumina_variants.tsv BC02.trimmed.sorted.bam refs/ZIKV_REF.fasta > BC02.variants.freqs.txt
python scripts/freqs.py --variants illumina_variants.tsv BC03.trimmed.sorted.bam refs/ZIKV_REF.fasta > BC03.variants.freqs.txt
