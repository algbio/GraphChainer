import sys
# python splits.py ../data/LRC/LRC.gfa xx.fastq

def LoadGfaGraph(filename):
    VL, E = {}, {}
    for line in open(filename).readlines():
        if line[0] == 'S':
            # S	92533	A
            i, s = line[1:].strip().split()
            VL[int(i)] = s
        elif line[0] == 'L':
            # L	104890	+	104892	+	0M
            li, lr, ri, rr, ov = line[1:].strip().split()
            li, ri = int(li), int(ri)
            if li not in E:
                E[li] = []
            E[li].append(ri)
    return VL, E

graph_filename = sys.argv[-2]
VL, E = LoadGfaGraph(graph_filename)

filename = sys.argv[-1]
lines = open(filename).readlines()
name = '_'.join(filename.split('.')[:-1])
folder = name + "_splits"
import os
import os.path
os.makedirs(folder, exist_ok=True)
count = 0
block_size = 1

names = {}

for i, line in enumerate(lines):
    if line[0] == '@':
        count += 1
        # if count % block_size == 1:
        out = open(folder + "/%05d.fastq" % count, 'w')
        for j in range(4):
            out.write(lines[i + j])
        # if count % block_size == 0:
        out.close()
        name = lines[i].split()[0][1:]
        names[count] = name
# if count % block_size != 0:
#     out.close()


# count = 10

import json

import subprocess
def exec(cmd):
    print(cmd)
    # process = subprocess.Popen(['time', cmd])
    # process.wait()
    # print(process.returncode)
    os.system(cmd)

sss = {}
header = ''
for i, line in enumerate(open('summary_chainpath.csv').readlines()):
    if i == 0:
        header = line.strip()
    else:
        name = line.split(',')[0]
        sss[name] = line.strip()


chainpaths = {}
lines = open('chainpaths.fasta').readlines()
for i, line in enumerate(lines):
    if line[0] == '>':
        name = line.strip()[1:]
        chainpaths[name] = lines[i + 1].strip()


out = open('vg_long_summary.csv', 'w')
out.write('name,length,vg_seq_len,vg_score,vg_time,vg_mm2_score,chainpath_mm2_score,' + header + '\n')
out.flush()


for i in range(1, count + 1):
    name = names[i]
    if name == '37ce64fb-9afb-3a8c-87bf-c3e0d0143e8c':
        break
    continue
    print('now', i, '/', count, name)
    fq_in = folder + "/%05d.fastq" % i
    gam_out = folder + "/%05d.gam" % i
    log = folder + "/%05d.log.txt" % i
    if not os.path.exists(gam_out):
        exec('(time vg map -M 1 -f ' + fq_in + ' -x x.xg -g x.gcsa > ' + gam_out + ") 2>" + log)
    json_out = folder + "/%05d.json" % i
    if not os.path.exists(json_out):
        exec('(time vg view -a ' + gam_out + ' -j > ' + json_out + ") 2>>" + log)
    fasta_out = folder + "/%05d.fasta" % i
    if not os.path.exists(fasta_out):
        exec('(time python getpaths.py ../data/LRC/LRC.gfa ' + json_out + ' ' + fasta_out + ") 2>>" + log)
    vg_minimap2_sam_out = folder + "/%05d.mm2.sam" % i
    # if not os.path.exists(vg_minimap2_sam_out):
    exec('(time minimap2 ' + fasta_out + ' ' + fq_in + ' -a -o ' + vg_minimap2_sam_out + ' 2> ' + folder + "/%05d.log.mm2.txt" % i + ") 2>>" + log)
    
    cp_fasta_out = folder + "/%05d.chainpath.fasta" % i
    # if not os.path.exists(cp_fasta_out):
    oo = open(cp_fasta_out, 'w')
    oo.write('>' + name + '\n')
    oo.write(chainpaths[name] + '\n')    
    cp_minimap2_sam_out = folder + "/%05d.cp.mm2.sam" % i
    # if not os.path.exists(cp_minimap2_sam_out):
    exec('(time minimap2 ' + cp_fasta_out + ' ' + fq_in + ' -a -o ' + cp_minimap2_sam_out + ' 2> ' + folder + "/%05d.log.mm2.cp.txt" % i + ") 2>>" + log)
    
    chainpath_mm2_score = 0
    for line in open(cp_minimap2_sam_out).readlines():
        for s in line.strip().split():
            if s.startswith('AS:i:'):
                chainpath_mm2_score = int(s.split(':')[-1])


    vg_mm2_score = 0
    for line in open(vg_minimap2_sam_out).readlines():
        for s in line.strip().split():
            if s.startswith('AS:i:'):
                vg_mm2_score = int(s.split(':')[-1])

    vg_time = ''
    for line in open(log).readlines():
        if 'user' in line:
            vg_time = line.split()[0][:-4]
            break
    

    line = open(json_out).readlines()[0]
    d = json.loads(line.strip())
    ids = []
    name = d['name']
    vg_score = 0
    if 'path' in d:
        path = d['path']['mapping']
        is_reverse = False
        if 'refpos' in d and 'is_reverse' in d['refpos']:
            is_reverse = d['refpos']['is_reverse']
        for e in path:
            id = e['position']['node_id']
            if 'is_reverse' in e['position']:
                if len(VL[int(id)]) > 1:
                    is_reverse = e['position']['is_reverse']
            if len(ids) == 0 or id != ids[-1]:
                ids.append(id)
        if 'score' in d:
            vg_score = d['score']
    ps = list(map(int, ids))
    seq = ''.join(VL[i] for i in ps)
    out.write(name + ',' + str(len(d['sequence'])) + ',' + str(len(seq)) + ',' 
        + str(vg_score) + ',' + str(vg_time) + ',' + str(vg_mm2_score) + ',' + str(chainpath_mm2_score) + ',' + sss[name] + '\n')
    out.flush()

out.close()