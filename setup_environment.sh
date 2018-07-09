conda create -n zika-isnv
source activate zika-isnv
conda install porechop bwa samtools pysam pyvcf
git clone git@github.com:artic-network/fieldbioinformatics.git && cd fieldbioinformatics && python setup.py install && cd ..
wget http://nanopore.s3.climb.ac.uk/Zika_iSNV_R9.4_albacore2.3.1.tgz
wget http://nanopore.s3.climb.ac.uk/Zika_iSNV_R9.5_albacore2.3.1.tgz
tar xvfz Zika_iSNV_R9.4_albacore2.3.1.tgz

