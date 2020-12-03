set -v

# time make -j 8 genread
# time ./genread -f ../data/LRC/LRC.gfa -o ref.fasta

# time badread simulate --reference ./ref.fasta --quantity 12x --length 15000,10000 \
#     --error_model pacbio --junk_reads 0 --random_reads 0 --chimeras 0 > badreads.fastq 2>badreads.log.txt

# vg index -x x.xg -g x.gcsa -k 16 x.vg
# python norev.py ref.fasta badreads.fastq badreads_norev.fastq
# python extract.py badreads_norev.fastq badreads.short.fasta
# rm badreads_norev.fastq.fai


# python -i vg_align_all.py badreads_norev.fastq



time python splits.py badreads_norev.fastq
time vg map -M 1 -f badreads_norev.fastq -x x.xg -g x.gcsa > vg_aln_long.gam
# time vg surject -x x.xg -b vg_aln_long.gam > vg_aln_long.bam

# time samtools view vg_aln_long.bam -o vg_aln_long.sam
time vg view -a badreads_norev_splits/vg_aln_long_00001.gam -j > badreads_norev_splits/vg_aln_long_00001.json
time python getpaths.py ../data/LRC/LRC.gfa badreads_norev_splits/vg_aln_long_00001.json badreads_norev_splits/vg_aln_long_00001.fasta

# time make -j 8 simple
# time ./simple -f ../data/LRC/LRC.gfa -a anchors.txt -o chains.txt

# time make -j 8 evaluate
# time ./evaluate -f ../data/LRC/LRC.gfa -p path.bin -c chains.txt -a anchors.txt -r badreads_norev.fastq -o summary.csv

