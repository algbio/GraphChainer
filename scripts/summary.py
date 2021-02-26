import sys
import os
import os.path
import json
import edlib
#            0              1                     2                  3    4
# python summary.py ../data/LRC/LRC.gfa ./data/badreads.fastq ga outdir

graph = sys.argv[1]
reads = sys.argv[2]
id = sys.argv[3]

outdir = "out"
if 4 < len(sys.argv):
	outdir = sys.argv[4]
if not outdir.endswith(id + '/'):
	outdir = outdir + '/' + id + '/'

def LoadGfaGraph(filename):
	VL, E = {}, {}
	for line in open(filename).readlines():
		if line[0] == 'S':
			# S 92533   A
			i, s = line[1:].strip().split()
			VL[int(i)] = s
		elif line[0] == 'L':
			# L 104890  +   104892  +   0M
			li, lr, ri, rr, ov = line[1:].strip().split()
			li, ri = int(li), int(ri)
			if li not in E:
				E[li] = []
			E[li].append(ri)
	return VL, E

VL, E = LoadGfaGraph(graph)

ed_global = lambda s1, s2 : edlib.align(s1, s2, mode='NW')['editDistance']
ed_local = lambda s1, s2 : edlib.align(s1, s2, mode='HW')['editDistance']
list2idx = lambda a: { a[i] : i for i in range(len(a)) }
revc = lambda s: ''.join({"A":"T","T":"A","C":"G","G":"C"}[c] for c in s[::-1])

import vg_pb2
import gzip
from google.protobuf.internal.encoder import _VarintBytes
from google.protobuf.internal.decoder import _DecodeVarint32, _VarintDecoder
def _VarintDecoder(mask):
	local_ord = ord
	def DecodeVarint(buffer, pos):
		result = 0
		shift = 0
		while 1:
			b = local_ord(buffer[pos])
			result |= ((b & 0x7f) << shift)
			pos += 1
			if not (b & 0x80):
				result &= mask
				return (result, pos)
			shift += 7
			# if shift >= 64:
			# 	raise _DecodeError('Too many bytes when decoding varint.')
	return DecodeVarint
_DecodeVarint64 = _VarintDecoder((1 << 64) - 1)
def read_alignments(gam_filename):
	with open(gam_filename, 'rb') as f:
		buf = gzip.GzipFile(fileobj=f).read()
		n = 0
		while n < len(buf):
			an, n = _DecodeVarint32(buf, n)
			for i in range(an):
				msg_len, n = _DecodeVarint32(buf, n)
				msg_buf = buf[n:n+msg_len]
				n += msg_len
				aln = vg_pb2.Alignment()
				aln.ParseFromString(msg_buf)
				yield aln

def parse_alignment(aln):
	# bps = sum(len(VL[x.position.node_id]) for x in a.path.mapping)
	name = aln.name.split()[0]
	seq = ''

	rev_cnt = 0
	for x in aln.path.mapping:
		idx = x.position.node_id
		ll = VL[idx]
		if x.position.is_reverse:
			rev_cnt += 1
			seq += revc(ll)
		else:
			seq += ll
	return {'name':name, 'seq':seq, 'path_cnt':len(aln.path.mapping), 'revcnt':rev_cnt, 'path_bps':len(seq)}

def parse_gam(filename):
	ret = {}
	for aln in read_alignments(filename):
		a = parse_alignment(aln)
		ret[a['name']] = a
	return ret

seqs_long = parse_gam(f'{outdir}{id}_long.gam')
seqs_clcs = parse_gam(f'{outdir}{id}_clc.gam')

def read_fastq(fastq_filename):
	reads_lines = open(fastq_filename).readlines()
	for i, line in enumerate(reads_lines):
		if line[0] == '@':
			info = reads_lines[i].strip()
			# name = info.split()[0][1:]
			seq = reads_lines[i + 1].strip()
			yield (info, seq)

seqs_read = {info.split()[0][1:] : (seq, info) for info, seq in read_fastq(f'{reads}')}

class CSV:
	def __init__(self):
		self.h = []
		self.hidx = {}
		self.r = []
		self.ridx = {}
		self.data = []
	def add_headers(self, headers):
		for h in headers:
			if h not in self.hidx:
				self.h.append(h)
				self.hidx[h] = len(self.h) - 1
	def get_hids(self, headers):
		return [self.hidx[x] for x in headers]
	def add(self, row, hids = []):
		if len(hids) == 0:
			self.data.append(row[:])
		else:
			tmp = [''] * len(self.h)
			for i in range(len(hids)):
				tmp[hids[i]] = row[i]
			self.data.append(tmp[:])
	def save(self, filename):
		fout = open(filename, 'w')
		fout.write(','.join(self.h) + '\n')
		for d in self.data:
			if len(d) < len(self.h):
				d += [''] * (len(self.h) - len(d))
			fout.write(','.join(d) + '\n')
		fout.close()

csv = CSV()
csv.add_headers(['name', 'length', 'br_id_rate']) #0,1,2
csv.add_headers(['long_pathcnt', 'long_path_bps', 'long_revcnt']) #3,4,5
csv.add_headers(['clcs_pathcnt', 'clcs_path_bps', 'clcs_revcnt']) #6,7,8
csv.add_headers(['long_align_rate']) #9
csv.add_headers([
'global_ed_read_long',  #10
# 'global_ed_long_true',
'global_ed_read_clcs',  #11
# 'global_ed_clcs_true',
# 'local_ed_read_long',
# 'local_ed_long_read',
# 'local_ed_true_long',
# 'local_ed_long_true',
# 'local_ed_read_clcs',
# 'local_ed_clcs_read',
# 'local_ed_true_clcs',
# 'local_ed_clcs_true',
# 'global_ed_read_clcs',
# 'global_ed_long_clcs',
# 'global_ed_read_true',
# 'global_ed_clcs_true',
])

reads_cnt = 0
for name in seqs_read:
	reads_cnt += 1
	if reads_cnt % (len(seqs_read) // 5 + 1) == 0:
		print(reads_cnt, '/', len(seqs_read))
	seq, info = seqs_read[name]
	row = [''] * len(csv.h)
	row[0] = name
	for t in info.split():
		if t.startswith('length='):
			row[1] = t.split('=')[-1]
	row[2] = str('%.3f'%(float(info.split()[-1].split('=')[-1][:-1]) / 100))

	long_seq = ''
	if name in seqs_long:
		a = seqs_long[name]
		long_seq = a['seq']
		row[3] = str(a['path_cnt'])
		row[4] = str(a['path_bps'])
		row[5] = str(a['revcnt'])
		row[10] = str(ed_global(seq, long_seq))
	row[9] = str(len(long_seq) / len(seq))
	clcs_seq = ''
	if name in seqs_clcs:
		a = seqs_clcs[name]
		clcs_seq = a['seq']
		row[6] = str(a['path_cnt'])
		row[7] = str(a['path_bps'])
		row[8] = str(a['revcnt'])
		row[11] = str(ed_global(seq, clcs_seq))
	csv.add(row)

csv.save(f'{outdir}{id}_summary.csv')
