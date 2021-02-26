# output cmd to screen before execution
set -x

# parameters
G=1000
L=150
S=33

# identifiers for this experiment
idx=LRC_15x

# data files
Graph=/mnt/c/Code/Summer/GCimplements/data/LRC/LRC.gfa
Reads=/mnt/d/summer/data/$idx.fastq

# path to binary
Bin=../bin/GraphAligner

# check whether dataset exists
if [ ! -f "$Reads" ]; then
    # sample a longest path on the graph as referrence sequence
    Ref=/mnt/d/summer/data/$idx.fasta
    time $Bin --generate-path -g $Graph -f $Ref -x vg -a /mnt/d/summer/data/$idx.gam
    # idx.txt has the node indices of the generated path
    time mv /mnt/d/summer/data/$idx.gam /mnt/d/summer/data/$idx.txt
    # simulate a PacBio long read dataset of 15x coverage using `badread` (commit 9e030e84849281e7dc92f0c9767b601c4dc9701e from https://github.com/rrwick/Badread.git)
    time badread simulate --reference $Ref --quantity 15x --length 15000,10000 --error_model pacbio --junk_reads 0 --random_reads 0 --chimeras 0 > $Reads 2>/mnt/d/summer/data/$idx\_br.log
fi

# experiment output folder
Out=out_$idx/ga
mkdir -p $Out/
Log=$Out/log.txt

echo "params=" $G $L $S > $Log

# use /usr/bin/time to measure "Max memory(kb)", "User time(s)", "System time(s)"
TimeCmd="/usr/bin/time -o $Log -a -f '%Mkb,%Us,%Ss'"

# align directly by original GraphAligner
$TimeCmd $Bin -x vg -f $Reads -g $Graph -a $Out/ga_long.gam 1>>$Log 2>>$Log

# align by colinear chaining
$TimeCmd $Bin --colinear-chaining -x vg -f $Reads -g $Graph -a $Out/ga_clc.gam --colinear-gap $G --colinear-split-len $L colinear-split-gap $S  1>>$Log 2>>$Log

# generate summary.csv
$TimeCmd python summary.py $Graph $Reads ga out_$idx  1>>$Log 2>>$Log


