import sys
# python splits.py xx.fastq
filename = sys.argv[-1]
lines = open(filename).readlines()
name = '_'.join(filename.split('.')[:-1])

import os
os.makedirs(name + "_splits", exist_ok=True)
count = 0
block_size = 1

for i, line in enumerate(lines):
    if line[0] == '@':
        count += 1
        if count % block_size == 1:
            out = open(name + "_splits/%05d.fastq" % count, 'w')
        for j in range(4):
            out.write(lines[i + j])
        if count % block_size == 0:
            out.close()
if count % block_size != 0:
    out.close()
