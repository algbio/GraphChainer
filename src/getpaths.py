import json
import sys
# python getpaths.py ../data/LRC/LRC.gfa vg_aln_long.json vg_aln_paths.fasta
graph_filename = sys.argv[-3]

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

VL, E = LoadGfaGraph(graph_filename)

def is_path(VL, E, p):
    for i in range(len(p)):
        if p[i] not in VL:
            print(p[i], 'not found')
            return False
        if i > 0 and p[i] not in E[p[i - 1]]:
            print(p[i-1], 'not connected to', p[i])
            return False
    return True

json_filename = sys.argv[-2]
out = open(sys.argv[-1], 'w')
lines = open(json_filename).readlines()

names = {}

for i, line in enumerate(lines):
    d = json.loads(line.strip())
    name = d['name']
    if name not in names:
        names[name] = 0
    names[name] += 1
    # ns = name.split('_')
    if 'path' not in d:
        continue
    # name, l, r = ns[0], ns[-2], ns[-1]
    path = d['path']['mapping']
    # out.write('>' + name + ' ' + l + ' ' + r + '\n')
    is_reverse = False
    if 'refpos' in d and 'is_reverse' in d['refpos']:
        is_reverse = d['refpos']['is_reverse']
    ids = []
    for e in path:
        id = e['position']['node_id']
        if 'is_reverse' in e['position']:
            if len(VL[int(id)]) > 1:
                is_reverse = e['position']['is_reverse']
                # if is_reverse:
                #     print(len(VL[int(id)]), e['position'])
        # offset = e['position']['offset']
        # l = e[]
        # if is_reverse:
        #     id = '-' + id
        if len(ids) == 0 or id != ids[-1]:
            ids.append(id)
        # print(ids, id)
    # break
    is_reverse = False
    if not is_reverse:
        ps = list(map(int, ids))
        # if not is_path(VL, E, ps):
        #     # print("? not a path ", name, names[name])
        #     continue
        seq = ''.join(VL[i] for i in ps)
        out.write('>' + name + "_%d"%names[name] + '\n')
        out.write(seq + '\n')
        # out.write(' '.join(ids) + '\n')
        # out.flush()
    print(name, names[name], is_reverse, len(seq), len(d['sequence']), 'score =', d['score'])
    # if len(seq) < len(d['sequence']) // 3:
    #     break
    
out.close()
