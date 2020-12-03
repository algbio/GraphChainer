import json
import sys
# python norev.py ref.fasta i.fastq out.fastq
filename = sys.argv[-2]
lines = open(filename).readlines()
outs = []

LEN = sum(len(l.strip()) for l in open(sys.argv[-3]).readlines()[1:])

i = 0
cp = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C' }
while i < len(lines):
    if lines[i][0] == '@':
        # @1854f7b0-ed60-7062-55d7-f62d4c68cc45 tmp,+strand,163965-184262 length=19947 error-free_length=20175 read_identity=93.93%
        if ',-strand,' in lines[i]:
            pass
            # ss = lines[i].split(' ')
            # tt = ss[1].split(',')
            # l, r = map(int, tt[2].split('-'))
            # # tt[1] = '+strand'
            # l, r = LEN - r + 1, LEN - l + 1
            # tt[2] = str(l) + '-' + str(r)
            # ss[1] = ','.join(tt)
            # lines[i] = ' '.join(ss)
            # lines[i + 1] = ''.join(cp[c] for c in lines[i + 1].strip())[::-1]
            # lines[i + 3] = lines[i + 3].strip()[::-1]
            # for j in range(4):
            #     outs.append(lines[i + j].strip() + '\n')
        else:
            for j in range(4):
                outs.append(lines[i + j].strip() + '\n')
        i += 4
    else:
        i += 1

open(sys.argv[-1], 'w').writelines(outs)