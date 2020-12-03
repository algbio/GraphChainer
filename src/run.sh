set -v

# time make -j 8 genread
# time ./genread -f ../data/LRC/LRC.gfa -o ref.fasta

# time badread simulate --reference ./ref.fasta --quantity 12x --length 15000,10000 \
#     --error_model pacbio --junk_reads 0 --random_reads 0 --chimeras 0 > badreads.fastq 2>badreads.log.txt

# vg index -x x.xg -g x.gcsa -k 16 x.vg
python norev.py ref.fasta badreads.fastq badreads_norev.fastq
python extract.py badreads_norev.fastq badreads.short.fasta
rm badreads.short.fasta.fai
time vg map -M 10 -F badreads.short.fasta -x x.xg -g x.gcsa > aln.gam
# time vg surject -x x.xg -b aln.gam > aln.bam

# time samtools view aln.bam -o aln.sam
time vg view -a aln.gam -j > aln.json
time python getanchors.py aln.json anchors.txt

time make -j 8 simple
time ./simple -f ../data/LRC/LRC.gfa -a anchors.txt -o chains.txt

time make -j 8 evaluate
time ./evaluate -f ../data/LRC/LRC.gfa -p path.bin -q chainpaths.fasta -c chains.txt -a anchors.txt -r badreads_norev.fastq -o summary.csv

