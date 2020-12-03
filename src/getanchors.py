import json
import sys
# python getchains.py aln.json chains.txt
filename = sys.argv[-2]
out = open(sys.argv[-1], 'w')
lines = open(filename).readlines()

names = set()

for i, line in enumerate(lines):
    d = json.loads(line.strip())
    name = d['name']
    # if name in names:
    #     print(name)
    names.add(name)
    ns = name.split('_')
    if 'path' not in d:
        continue
    name, l, r = ns[0], ns[-2], ns[-1]
    path = d['path']['mapping']
    out.write('>' + name + ' ' + l + ' ' + r + '\n')
    is_reverse = False
    if 'is_reverse' in d['refpos']:
        is_reverse = d['refpos']['is_reverse']
    ids = []
    for e in path:
        if 'is_reverse' in e['position']:
            is_reverse = e['position']['is_reverse']
        id = e['position']['node_id']
        # offset = e['position']['offset']
        # l = e[]
        if is_reverse:
            id = '-' + id
        ids.append(id)
    out.write(' '.join(ids) + '\n')
    out.flush()
    
out.close()
