cd ../../data/css/raw/pa-seq

for file in `ls *.sra`
do
	bsub ~/sra/bin64/fastq-dump -O fastq --gzip --split-3 $file
done
