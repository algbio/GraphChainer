
import sys
# python extract.py in.fastq out.filename
filename = sys.argv[-2]
# name = '.'.join(filename.split('.')[:-1])
out = open(sys.argv[-1], 'w')
lines = open(filename).readlines()
L = 150
offset = 36
for i, line in enumerate(lines):
    if line[0] == '@':
        seq = lines[i + 1].strip()
        N = (len(seq) - L + 1) // offset + 1
        for j in range(N):
            l = j*(len(seq) - L + 1)//N
            s = seq[l : l + L]
            out.write('>' + line[1:].split()[0] + '__short_%d_%d ' % (l, l + L) + ' '.join(line.split()[1:]) + '\n')
            out.write(s + '\n')
out.close()
