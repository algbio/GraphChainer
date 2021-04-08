

import os

N = 2
Bin = "./bin/GraphAligner"
Graphs = ["/mnt/c/Code/Summer/GCimplements/data/LRC/LRC.vg", "/mnt/c/Code/Summer/GCimplements/data/MHC/MHC1.vg"]
Data = "./data/"
Gams = "./gams/"
Logs = "./logs/"
force_redo = True
Threads = 4
# default badread is length ~ (mean=15000,std=10000)
# default pbsim clr is length ~ (mean=3000,std=2300)
# or more `real` length ~ (mean=15000,std=10000)
read_length_mean, read_length_std = 3000, 2300

def mkdir_safe(path):
    if not os.path.exists(path):
        os.system(f"mkdir -p {path}")
mkdir_safe(Data)
mkdir_safe(Gams)
mkdir_safe(Logs)



def time_cmd(cmd, log):
    return f"/usr/bin/time -o {log} -a -v {cmd}"
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
params.append((150,150,10000))
# [(100,31), (100,17), ]
# for L, S in [(150,150)]:
#     for G in [10000]:
#         params.append((L, S, G))

gen_log = f"{Logs}/gen.log.txt"
for Graph in Graphs:
    for idx in range(N):
        seed = idx
        id = Graph.split('/')[-1].split('.')[0] + '_' + str(idx)
        Reads_prefix = f"{Data}/{id}"
        Reads = f"{Reads_prefix}.fastq"
        if not os.path.exists(Reads):
            print(f"generating read set #{idx} to {id}")
            Ref = f"{Data}/{id}.fasta"
            run(f"{Bin} --generate-path --generate-path-seed {seed} -g {Graph} -f {Ref} -x vg -a {Data}/{id}.gam", gen_log)
            # id.path.txt has the node indices of the generated path
            run(f"mv {Data}/{id}.gam {Data}/{id}.path.txt", gen_log)
            # # simulate a PacBio long read dataset of 15x coverage using `badread` (commit 9e030e84849281e7dc92f0c9767b601c4dc9701e from https://github.com/rrwick/Badread.git)
            # run(f"badread simulate --seed {seed} --reference {Ref} --quantity 15x --length 15000,10000 --error_model pacbio --junk_reads 0 --random_reads 0 --chimeras 0 > {Reads} 2>{Data}/{id}_br.log.txt", gen_log, True, False)
            
            # simulate a PacBio long read dataset of 20x coverage using 'pbsim' (commit e014b1dd40e87a8799346a9835d70a4da3dc857c from https://github.com/pfaucon/PBSIM-PacBio-Simulator.git)
            prefix = 'xxxx'
            run(f"./bin/pbsim --data-type CLR --depth 20 --seed {seed} --model_qc ./bin/model_qc_clr {Ref} --prefix {prefix} --length-mean {read_length_mean} --length-sd {read_length_std}", gen_log, True, True)
            # rename the reads generated
            run(f"mv ./{prefix}_0001.fastq {Reads}")
            ReadsMaf = f"{Reads_prefix}.maf"
            run(f"mv ./{prefix}_0001.maf {ReadsMaf}")

            # run(f"badread simulate --seed {seed} --reference {Ref} --quantity 15x --length 15000,10000 --error_model pacbio --junk_reads 0 --random_reads 0 --chimeras 0 > {Reads} 2>{Data}/{id}_br.log.txt", gen_log, True, False)
        
        long_log = f"{Logs}/{id}_long.log.txt"
        long_gam = f"{Gams}/{id}_long.gam"
        if force_redo or not os.path.exists(long_gam):
            run(f"{Bin} -t {Threads} -x vg -f {Reads} -g {Graph} -a {long_gam}", long_log)
        
        aln_clc_log = f"{Logs}/{id}_long.log.txt"
        for L, S, G in params:
            clc_gam = f"{Gams}/{id}_clc_{L}_{S}_{G}.gam"
            clc_log = f"{Gams}/{id}_clc_{L}_{S}_{G}.log.txt"
            if force_redo or not os.path.exists(clc_gam):
                run(f"{Bin} -t {Threads} --colinear-chaining -x vg -f {Reads} -g {Graph} -a {clc_gam} --colinear-gap {G} --colinear-split-len {L} --colinear-split-gap {S} --short-verbose", clc_log)



