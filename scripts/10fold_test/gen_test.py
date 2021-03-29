

import os

N = 1
Bin = "./bin/GraphAligner"
Graph = "/mnt/c/Code/Summer/GCimplements/data/LRC/LRC.gfa"
Data = "/mnt/d/summer/data/"
Gams = "/mnt/d/summer/gams/"

if not os.path.exists(Gams):
    os.system(f"mkdir -p {Gams}")

def time_cmd(cmd, log):
    return f"/usr/bin/time -o {log} -a -f 'TIME:%Mkb,%Us,%Ss' {cmd}"
def log_cmd(cmd, log):
    return f"{cmd} 1>>{log} 2>>{log}"
def run(cmd, log = "", t = True, l = True, v = True):
    if v:
        print(cmd)
    if log != "":
        open(log, 'a').write(cmd + '\n')
        if t:
            cmd = time_cmd(cmd, log)
        if l:
            cmd = log_cmd(cmd, log)
    try:
        os.system(cmd)
    except KeyboardInterrupt:
        exit(1)

params = []
# [(100,31), (100,17), ]
for L, S in [(150, 31), (150, 71), (150,150), (300, 71), (300, 300)]:
    for G in [1000, 10000]:
        params.append((L, S, G))

gen_log = "./gen.log.txt"
for idx in range(N):
    seed = idx
    id = Graph.split('/')[-1].split('.')[0] + '_' + str(idx)
    Reads = f"{Data}/{id}.fastq"
    if not os.path.exists(Reads):
        print(f"generating read set #{idx} to {id}")
        Ref = f"{Data}/{id}.fasta"
        run(f"{Bin} --generate-path --generate-path-seed {seed} -g {Graph} -f {Ref} -x vg -a {Data}/{id}.gam", gen_log)
        # id.path.txt has the node indices of the generated path
        run(f"mv {Data}/{id}.gam {Data}/{id}.path.txt", gen_log)
        # simulate a PacBio long read dataset of 15x coverage using `badread` (commit 9e030e84849281e7dc92f0c9767b601c4dc9701e from https://github.com/rrwick/Badread.git)
        run(f"badread simulate --seed {seed} --reference {Ref} --quantity 15x --length 15000,10000 --error_model pacbio --junk_reads 0 --random_reads 0 --chimeras 0 > {Reads} 2>{Data}/{id}_br.log.txt", gen_log, True, False)
    
    aln_log = f"./{id}.log.txt"
    long_gam = f"{Gams}/{id}_long.gam"
    if not os.path.exists(long_gam):
        run(f"{Bin} -x vg -f {Reads} -g {Graph} -a {long_gam}", aln_log)
    
    for L, S, G in params:
        clc_gam = f"{Gams}/{id}_clc_{L}_{S}_{G}.gam"
        if not os.path.exists(clc_gam):
            run(f"{Bin} --colinear-chaining -x vg -f {Reads} -g {Graph} -a {clc_gam} --colinear-gap {G} --colinear-split-len {L} --colinear-split-gap {S} ", aln_log)







# # output cmd to screen before execution
# set -x

# # parameters
# G=1000
# L=150
# S=33

# # identifiers for this experiment
# idx=LRC_15x

# # data files
# Graph=/mnt/c/Code/Summer/GCimplements/data/LRC/LRC.gfa
# Reads=/mnt/d/summer/data/$idx.fastq

# # path to binary
# Bin=./bin/GraphAligner

# # check whether dataset exists
# if [ ! -f "$Reads" ]; then
#     # sample a longest path on the graph as referrence sequence
#     Ref=/mnt/d/summer/data/$idx.fasta
#     time $Bin --generate-path -g $Graph -f $Ref -x vg -a /mnt/d/summer/data/$idx.gam
#     # idx.txt has the node indices of the generated path
#     time mv /mnt/d/summer/data/$idx.gam /mnt/d/summer/data/$idx.txt
#     # simulate a PacBio long read dataset of 15x coverage using `badread` (commit 9e030e84849281e7dc92f0c9767b601c4dc9701e from https://github.com/rrwick/Badread.git)
#     time badread simulate --reference $Ref --quantity 15x --length 15000,10000 --error_model pacbio --junk_reads 0 --random_reads 0 --chimeras 0 > $Reads 2>/mnt/d/summer/data/$idx\_br.log
# fi

# # experiment output folder
# OutId=ga2
# Out=out_$idx/$OutId
# mkdir -p $Out/
# Log=$Out/log.txt

# echo "params=" $G $L $S > $Log

# # use /usr/bin/time to measure "Max memory(kb)", "User time(s)", "System time(s)"
# TimeCmd="/usr/bin/time -o $Log -a -f '%Mkb,%Us,%Ss'"

# # align directly by original GraphAligner
# $TimeCmd $Bin -x vg -f $Reads -g $Graph -a $Out/ga_long.gam 1>>$Log 2>>$Log

# # align by colinear chaining
# $TimeCmd $Bin --colinear-chaining -x vg -f $Reads -g $Graph -a $Out/ga_clc.gam --colinear-gap $G --colinear-split-len $L colinear-split-gap $S  1>>$Log 2>>$Log

# # generate summary.csv
# $TimeCmd python summary.py $Graph $Reads $OutId out_$idx  1>>$Log 2>>$Log


