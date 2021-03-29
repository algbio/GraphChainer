import sys
import os
import os.path
import json
import edlib
import time
import statistics
import numpy as np


Graph = "/mnt/c/Code/Summer/GCimplements/data/LRC/LRC.gfa"
Data = "/mnt/d/summer/data/"
Gams = "/mnt/d/summer/gams/"
Csvs = "/mnt/d/summer/csvs/"
Soss = "/mnt/d/summer/csvs/sos_cache/"
verbose = True
force_redo_output = True

if not os.path.exists(Csvs):
    os.system(f"mkdir -p {Csvs}")
if not os.path.exists(Soss):
    os.system(f"mkdir -p {Soss}")

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

VL, E = LoadGfaGraph(Graph)


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
	# return {'name':name, 'seq':seq, 'path_cnt':len(aln.path.mapping), 'revcnt':rev_cnt, 'path_bps':len(seq)}
	return (name, seq, len(aln.path.mapping), rev_cnt, len(seq))

def parse_gam(filename):
	ret = {}
	for aln in read_alignments(filename):
		a = parse_alignment(aln)
		ret[a[0]] = a
		# ret[a['name']] = a
	return ret

gams_all = {}
def lazy_load_gam(name):
	if name not in gams_all:
		gams_all[name] = parse_gam(f"{Gams}/{name}.gam")
	return gams_all[name]

def read_fastq(fastq_filename):
	reads_lines = open(fastq_filename).readlines()
	for i, line in enumerate(reads_lines):
		if line[0] == '@':
			info = reads_lines[i].strip()
			# name = info.split()[0][1:]
			seq = reads_lines[i + 1].strip()
			yield (info, seq)

datasets_all = {}
def lazy_load_fastq(name):
	if name not in datasets_all:
		reads = {info.split()[0][1:] : (seq, info) for info, seq in read_fastq(f"{Data}/{name}.fastq")}
		ref = open(f"{Data}/{name}.fasta").readlines()[-1].strip()
		datasets_all[name] = (reads, ref)
	return datasets_all[name]

log_time_all = {}
def lazy_load_time(name):
	if name not in log_time_all:
		# LRC_0_clc_150_31_1000.gam
		id = '_'.join(name.split('_')[:2])
		# graph, idx, _, L, S, G = name.split('_')
		# id = graph + '_' + str(idx)
		lines = open(f"./{id}.log.txt").readlines()
		now = None
		for i, line in enumerate(lines):
			if line.startswith("Co-linear chaining"):
				if line.startswith("Co-linear chaining on splits="):
					# Co-linear chaining on splits=(150,31,1000)
					l, s, g = line.split('(')[1].split(')')[0].split(',')
					now = f'{id}_clc_{l}_{s}_{g}'
				elif line.startswith("Co-linear chaining off"):
					now = f"{id}_long"
			elif line.startswith("TIME:"):
				# TIME:254064kb,33.45s,0.45s
				pm, ut, st = line.strip().split(':')[1].split(',')
				pm, ut, st = pm[:-2], ut[:-1], st[:-1]
				pm = '%.2f' % (float(pm) / 1024)
				if not (now is None):
					log_time_all[now] = (pm, ut, st)
				now = None
	return log_time_all[name]


eds_all = {}
def lazy_load_eds(name):
	#name is LRC_1_long or LRC_0_clc_300_71_10000
	if name not in eds_all:
		eds = {}
		cache = f"{Soss}/{name}_ed.csv"
		if os.path.exists(cache):
			if verbose:
				print('load eds', name)
			for line in open(cache).readlines()[1:]:
				ss = line.strip().split(',')
				id = ss[0]
				row = list(map(int, ss[2:]))
				eds[id] = row

		else:
			if verbose:
				print('compute eds', name)
			gam = lazy_load_gam(name)
			reads_name = "_".join(name.split('_')[:2])
			reads, ref_seq = lazy_load_fastq(reads_name)
			cf = open(cache, 'w')
			cf.write('id,len,aln_len,ed_read,ed_true\n')
			for id in reads:
				read_seq, info = reads[id]
				true_seq = ''
				s, t = map(int, info.split()[1].split(',')[-1].split('-'))
				if '-strand' in info.split()[1].split(','):
					true_seq = revc(ref_seq[len(ref_seq) - t : len(ref_seq) - s + 1])
				else:
					true_seq = ref_seq[s : t + 1]

				long_seq = ''
				if id in gam:
					a = gam[id]
					long_seq = a[1] #a['seq']
				row = []
				row.append(len(long_seq))
				row.append(ed_global(read_seq, long_seq))
				row.append(ed_global(true_seq, long_seq))

				eds[id] = row
				cf.write(','.join([id] + list(map(str, [len(read_seq)] + row))) + '\n')
			cf.close()

		eds_all[name] = eds
	return eds_all[name]

if not os.path.exists(f"{Csvs}/sos.csv") or force_redo_output:
	f = open(f"{Csvs}/sos.csv", 'w')
	header = 'id,grpah,L,S,G, selected,comment'
	header += ', peak_memory(mb),user_time(s),system_time(s)'
	header += ', rough_better_cnt'
	header += ', avg_len(aln)/len(read),std,mid,l0.2,r0.2'
	header += ', avg_ed(aln_read)/len(read),std,mid,l0.2,r0.2'
	header += ', avg_ed(aln_true)/len(read),std,mid,l0.2,r0.2'
	f.write(header + '\n')
	f.close()
sos_out = open(f"{Csvs}/sos.csv", 'a')

def asmlr(vals):
	n = len(vals)
	a = sum(vals) / n
	s = statistics.stdev(vals)
	t = sorted(vals)
	m = t[n//2]
	l = t[int(n * 0.2)]
	r = t[int(n * 0.8)]
	return a, s, m, l, r
done = set()
def analyze(f):
	if verbose:
		print('analyze', f)
	# LRC_0_clc_150_31_1000.gam
	name =  f.split('/')[-1].split('.')[0]
	reads_name = "_".join(name.split('_')[:2])
	reads, ref_seq = lazy_load_fastq(reads_name)
	graph, idx, _, L, S, G = name.split('_')
	long_name = f"{graph}_{idx}_long"
	long_eds = lazy_load_eds(long_name)
	clcs_eds = lazy_load_eds(name)
	ids_list = []
	ids_list.append(('all', [id for id in reads]))
	for e, l in [(0.3, 0.9), (0.2, 0.9)]:
		cmt = 'ed>%.1f|l<%.1f' % (e, l)
		ids_list.append((cmt, [id for id in reads if long_eds[id][1] / len(reads[id][0]) > e or long_eds[id][0] < len(reads[id][0]) * l]))
	
	def get_vals(eds, ids):
		leds, reds, teds = [], [], []
		rough_better_cnt = 0
		for id in ids:
			ll = len(reads[id][0])
			l, r, t = eds[id]
			leds.append(l / ll)
			reds.append(r / ll)
			teds.append(t / ll)
			if r < long_eds[id][1] - ll * 0.01:
				rough_better_cnt += 1
		return [rough_better_cnt] + list(['%.5f'%f for f in [*asmlr(leds), *asmlr(reds), *asmlr(teds)]])

	if long_name not in done:
		for comment, ids in ids_list:
			row = [long_name, graph, 0,0,0, len(ids),comment]
			row += list(lazy_load_time(long_name))
			row += get_vals(long_eds, ids)
			sos_out.write(','.join(map(str, row)) + '\n')
			sos_out.flush()
			if verbose:
				print(','.join(map(str, row)))
		done.add(long_name)

	for comment, ids in ids_list:
		row = [name, graph, L,S,G, len(ids),comment]
		row += list(lazy_load_time(name))
		row += get_vals(clcs_eds, ids)
		sos_out.write(','.join(map(str, row)) + '\n')
		sos_out.flush()
		if verbose:
			print(','.join(map(str, row)))

if True:
	todo = []
	for f in os.listdir(f"{Gams}"):
		if 'clc' not in f:
			continue
		if f not in done:
			todo.append((f, os.stat(f"{Gams}/{f}").st_size))
	time.sleep(1)
	for f, sz in todo:
		tsz = os.stat(f"{Gams}/{f}").st_size
		if tsz == sz:
			analyze(f)
			done.add(f)

sos_out.close()
# exit(0)
#### short sos
sos_lines = open(f"{Csvs}/sos.csv").readlines()
old_header = 'id,grpah,L,S,G, selected,comment'
old_header += ', peak_memory(mb),user_time(s),system_time(s)'
old_header += ', rough_better_cnt'
old_header += ', avg_len(aln)/len(read),avg_len(aln)/len(read)_std,avg_len(aln)/len(read)_mid,avg_len(aln)/len(read)_l0.2,avg_len(aln)/len(read)_r0.2'
old_header += ', avg_ed(aln_read)/len(read),avg_ed(aln_read)/len(read)_std,avg_ed(aln_read)/len(read)_mid,avg_ed(aln_read)/len(read)_l0.2,avg_ed(aln_read)/len(read)_r0.2'
old_header += ', avg_ed(aln_true)/len(read),avg_ed(aln_true)/len(read)_std,avg_ed(aln_true)/len(read)_mid,avg_ed(aln_true)/len(read)_l0.2,avg_ed(aln_true)/len(read)_r0.2'
old_header = old_header.replace(' ', '')
old_cols = old_header.split(',')
pars = {}
for line in sos_lines[1:]:
	ss = line.strip().split(',')
	name = ss[0]
	id = name.split('_')[0] + '_' + '_'.join(name.split('_')[2:])
	id += '_' + ss[old_cols.index('comment')]
	# LRC_0_clc_150_31_1000.gam
	if id not in pars:
		pars[id] = []
	pars[id].append(ss)

sosos_out = open(f"{Csvs}/sosos.csv", 'w')
same_cols = 'comment,grpah,L,S,G'.split(',')
calc_header = 'peak_memory(mb),user_time(s),system_time(s)'
calc_header += ',selected,rough_better_cnt'
calc_header += ',avg_len(aln)/len(read),avg_len(aln)/len(read)_std'
calc_header += ',avg_ed(aln_read)/len(read),avg_ed(aln_read)/len(read)_std'
calc_header += ',avg_ed(aln_true)/len(read),avg_ed(aln_true)/len(read)_std'
calc_header = calc_header.replace(' ', '')
calc_cols = calc_header.split(',')
new_header = 'id,trials' + ','.join(same_cols) + ','+ ','.join(['avg_'+s+',std_'+s for s in calc_cols])
sosos_out.write(new_header + '\n')
hs = new_header.strip().split(',')
for id in pars:
	vs = pars[id]
	n = len(vs)
	row = [id]
	for h in same_cols:
		row.append(vs[0][old_cols.index(h)])
	for h in calc_cols:
		idx = old_cols.index(h)
		vals = [float(v[idx]) for v in vs]
		a, s, m, l, r = asmlr(vals)
		row.append('%.5f'%a)
		row.append('%.5f'%s)
	sosos_out.write(','.join(row) + '\n')
sosos_out.close()

#### xxplots
import matplotlib.pyplot as plt
Pngs = "/mnt/d/summer/pngs/"
if not os.path.exists(Pngs):
    os.system(f"mkdir -p {Pngs}")
sepstrs = ['0.3', '0.1', '0.01']
seps = [float(s) for s in sepstrs]
plots_done = set()
stys = ['bo', 'rx', 'g^', 'pd']
mem = {}
for g in done:
	if '_long' in g:
		continue
	graph, idx, _, L, S, G = g.split('_')
	target = f'{L}_{S}_{G}'
	if target in plots_done:
		continue
	plots_done.add(target)
	allx, ally = [], []
	xs, ys = [[] for s in seps], [[] for s in seps]
	rs = []
	ws = []
	ls = []
	es = []
	r1s, r2s = [], []
	total_reads = 0
	all_eds = []
	for f in done:
		if target not in f:
			continue
		name =  f.split('/')[-1].split('.')[0]
		reads_name = "_".join(name.split('_')[:2])
		reads, ref_seq = lazy_load_fastq(reads_name)
		graph, idx, _, L, S, G = name.split('_')
		long_name = f"{graph}_{idx}_long"
		long_eds = lazy_load_eds(long_name)
		clcs_eds = lazy_load_eds(name)
		#         0          1          2           3             4       5           6              7              
		# [[ len(read),  len(aln),ed(aln,read),ed(aln,true),  len(aln),ed(aln,read),ed(aln,true),  errs ]]
		eds = np.ndarray((len(reads), 8)) 
		eds_idx = 0
		for id in reads:
			total_reads += 1
			ll = len(reads[id][0])
			l, r, t = clcs_eds[id]
			now = [ll] + clcs_eds[id] + long_eds[id]
			# leds.append(l / ll)
			# reds.append(r / ll)
			# teds.append(t / ll)
			# better = False
			# ls.append(long_eds[id][1] / ll)
			# rs.append(long_eds[id][1] / ll - r / ll)
			# r1s.append(long_eds[id][1] / ll)
			# r2s.append(r / ll)
			# ws.append(min(long_eds[id][1] / ll, r / ll))
			now.append(float( reads[id][1].split()[-1][:-1].split('=')[1] ) / 100.0)
			eds[eds_idx] = np.array(now)
			eds_idx += 1
			# for si in range(len(seps)):
			# 	if seps[si] < long_eds[id][1] / ll - r / ll:
			# 		xs[si].append(r / ll)
			# 		ys[si].append(l / ll)
			# 		better = True
			# 		break
			# if not better:
			# 	allx.append(r / ll)
			# 	ally.append(l / ll)
		all_eds.append(eds)
	
	eds = np.concatenate(all_eds)
	for i in range(1, 7):
		eds[:, i] /= eds[:, 0]

	# simple plots of curve of errors

	old_eds = eds[:, 5]
	new_eds = eds[:, 2]
	def plot_improvement(old_eds, new_eds):
		fig, ax1 = plt.subplots()
		ax1.set_xlabel('eds/len improvement')
		ax1.set_ylabel('number of reads')

		r_ips = old_eds - new_eds
		r_iids = np.argsort(r_ips)[::-1]
		maxx = r_ips.max()
		N = (r_ips > 0).sum()
		xs = list(range(1, N + 1))
		ax1.plot(r_ips[r_iids[:N]], xs, marker='o', markersize=5, label='improved')

		ax1.set_xlim(0, maxx)
		ax1.set_ylim(0, N)
		ax1.set_yticks(list(ax1.get_yticks()) + [N])
		ax1.legend(loc='upper left')

		ax2 = ax1.twinx()
		ax2.set_ylabel('exact eds/len')
		ax2.plot(r_ips[r_iids[:N]], old_eds[r_iids[:N]], marker='.', markersize=1, linewidth=1, alpha=0.2, color='orange', label='old')
		ax2.plot(r_ips[r_iids[:N]], old_eds[r_iids[:N]], marker='.', markersize=1, linewidth=0, alpha=0.7, color='orange')
		ax2.plot(r_ips[r_iids[:N]], new_eds[r_iids[:N]], marker='.', markersize=1, linewidth=1, alpha=0.2, color='green', label='clc')
		ax2.plot(r_ips[r_iids[:N]], new_eds[r_iids[:N]], marker='.', markersize=1, linewidth=0, alpha=0.7, color='green')
		ax2.legend(loc='upper right')
		return fig, ax1, ax2

	fig, ax1, ax2 = plot_improvement(eds[:, 5], eds[:, 2])
	plt.title(target.split('.')[0] + ' read edit distance improvement')
	plt.savefig(f'{Pngs}/{target}_all_rs_curve.png')
	plt.clf()

	fig, ax1, ax2 = plot_improvement(eds[:, 6], eds[:, 3])
	plt.title(target.split('.')[0] + ' true edit distance improvement')
	plt.savefig(f'{Pngs}/{target}_all_ts_curve.png')
	plt.clf()

	# if target == '150_150_10000.gam':
	# 	break

	# ax1.plot(m[0], m[1], '.', label='len_data')
	# ax1.plot(m[0], m[2], '.', label='len_content')
	# ax1.legend()

	# ax2 = ax1.twinx()
	# ax2.set_ylabel('ram')
	# ax2.plot(m[0], m[3], label='mem_start')
	# ax2.plot(m[0], m[4], label='mem_end')
	# ax2.legend()
	# plt.show()

	# ids = [(rs[i], r1s[i], r2s[i], i) for i in range(len(rs))]
	# ids.sort()
	# rs = sorted(rs)[::-1]
	# rs = [r for r in rs if r > 0]
	# N = len([r for r in rs if r > 0])
	# xs = [i for i in range(len(rs))]
	# plt.plot(rs, xs, marker='o', markersize=10, label='diff')
	# plt.plot([x[1] for x in ids if x[0] > 0], xs, label='old')
	# plt.plot([x[2] for x in ids if x[0] > 0], xs, label='new')
	# plt.ylabel('rs')
	# plt.xlabel('ids')
	
	# # better.png
	# ss = [len(x) for x in xs]
	# for i in range(1, len(ss)):
	# 	ss[i] += ss[i - 1]
	# plt.figure(figsize=(7,7))
	# for si in range(len(seps)):
	# 	plt.plot(xs[si], ys[si], color=stys[si%len(stys)][0], marker=stys[si%len(stys)][1], linestyle="None", label=f'>{sepstrs[si]}', alpha=0.2)
	# plt.plot(allx, ally, color='black', marker='.', linestyle="None", markersize=2, label=f'other', alpha=0.05)
	# plt.xlabel('ed_read_div_len')
	# plt.ylabel('len_aln_div_len')
	# plt.legend()
	# plt.xlim(0, 1.2)
	# plt.ylim(0, 1.2)
	# plt.savefig(f'{Pngs}/{target}_better_all.png')
	# plt.clf()
	
	# #separated error curve
	# # rs = sorted(rs)[::-1]
	# xs = [i / len(ls) for i in range(len(rs))]
	# ls = sorted([ls[i]-(1-es[i]) for i in range(len(ls))])[::-1]
	# ws = sorted([ws[i]-(1-es[i]) for i in range(len(ls))])[::-1]
	# plt.plot(ls, xs, label='GraphAligner')
	# plt.plot(ws, xs, label='GraphAligner+ColinearChaining')
	# xds, xxs = [], []
	# lsi, wsi = 0, 0
	# for x in sorted(list(set(ls + ws)))[::-1]:
	# 	while lsi + 1 < len(ws) and ls[lsi] >= x:
	# 		lsi += 1
	# 	while wsi + 1 < len(ws) and ws[wsi] >= x:
	# 		wsi += 1
	# 	xxs.append(x)
	# 	xds.append(xs[lsi] - xs[wsi])
	# plt.plot(xds, xxs, label='diff')
	# # plt.plot(xs, es, label='simulated error rate')
	# # plt.plot(xs, [e+0.1 for e in es], label='simulated error rate+0.1')
	# # plt.plot(xs, [e+0.2 for e in es], label='simulated error rate+0.2')
	# # plt.plot(xs, [e+0.3 for e in es], label='simulated error rate+0.3')
	# plt.ylabel('reads%%')
	# plt.xlabel('rate')
	# # plt.xticks()
	# plt.xlim(0, 0.5)
	# plt.ylim(0, 0.1)
	# # plt.yticks(list(plt.yticks()[0]) + [N])
	# plt.legend()
	# plt.savefig(f'{Pngs}/{target}_all_eder_curve.png')
	# plt.clf()


	# # compare with "simulated error rate"
	# # rs = [r for r in rs if r > 0]
	# # xs = [i for i in range(len(rs))]
	# # plt.plot(rs, xs, marker='o', markersize=10)
	# ls = sorted(ls)[::-1]
	# ws = sorted(ws)[::-1]
	# es = sorted([1 - e for e in es])[::-1]
	# plt.plot(xs, ls, label='GraphAligner')
	# plt.plot(xs, ws, label='GraphAligner+ColinearChaining')
	# plt.plot(xs, es, label='simulated error rate')
	# plt.plot(xs, [e+0.1 for e in es], label='simulated error rate+0.1')
	# plt.plot(xs, [e+0.2 for e in es], label='simulated error rate+0.2')
	# plt.plot(xs, [e+0.3 for e in es], label='simulated error rate+0.3')
	# plt.ylabel('rs')
	# plt.xlabel('ids')
	# # plt.xticks()
	# plt.xlim(0, 0.3)
	# N = len([r for r in rs if r > 0])
	# plt.ylim(0, 0.5)
	# # plt.yticks(list(plt.yticks()[0]) + [N])
	# plt.legend()
	# plt.savefig(f'{Pngs}/{target}_all_ed_curve.png')
	# plt.clf()


	# ids = [(rs[i], r1s[i], r2s[i], i) for i in range(len(rs))]
	# ids.sort()
	# rs = sorted(rs)[::-1]
	# rs = [r for r in rs if r > 0]
	# N = len([r for r in rs if r > 0])
	# xs = [i for i in range(len(rs))]
	# plt.plot(rs, xs, marker='o', markersize=10, label='diff')
	# plt.plot([x[1] for x in ids if x[0] > 0], xs, label='old')
	# plt.plot([x[2] for x in ids if x[0] > 0], xs, label='new')
	# plt.ylabel('rs')
	# plt.xlabel('ids')
	# plt.xlim(0, 1.0)
	# plt.ylim(0, N)
	# plt.yticks(list(plt.yticks()[0]) + [N])
	# plt.legend()
	# plt.savefig(f'{Pngs}/{target}_all_rs_curve.png')
	# plt.clf()
	
	# mem[target + '_rs'] = rs[:]
	# mem[target + '_ws'] = ws[:]
	# mem[target + '_ls'] = ls[:]
	# ff = lambda v: len([ x for x in v if x < 0.10 ])
	# print(target, ss, '/', total_reads, ff(ls), ff(ws))
	print(target, ss, '/', total_reads)




	# for i in range(1, len(ss)):
	# 	ss[i] += ss[i - 1]
	# print(target, ss, '/', total_reads)
	# # plt.figure(figsize=(20,20))
	# for si in range(len(seps)):
	# 	plt.scatter(xs[si], ys[si], stys[si%len(stys)], label=f'>{sepstrs[si]}', alpha=0.2)
	# plt.scatter(allx, ally, color='black', marker='.', markersize=1, label=f'other', alpha=0.1)
	# plt.xlabel('ed_read_div_len')
	# plt.ylabel('len_aln_div_len')
	# plt.legend()
	# # plt.xlim(0, 1.3)
	# # plt.ylim(0, 1.5)
	# plt.savefig(f'{Pngs}/{target}_better.png')
	# plt.clf()


# gen_log = "./gen.log.txt"
# for idx in range(N):
#     seed = idx
#     id = Graph.split('/')[-1].split('.')[0] + '_' + str(idx)
#     Reads = f"{Data}/{id}.fastq"
#     if not os.path.exists(Reads):
#         print(f"generating read set #{idx} to {id}")
#         Ref = f"{Data}/{id}.fasta"
#         run(f"{Bin} --generate-path --generate-path-seed {seed} -g {Graph} -f {Ref} -x vg -a {Data}/{id}.gam", gen_log)
#         # id.path.txt has the node indices of the generated path
#         run(f"mv {Data}/{id}.gam {Data}/{id}.path.txt", gen_log)
#         # simulate a PacBio long read dataset of 15x coverage using `badread` (commit 9e030e84849281e7dc92f0c9767b601c4dc9701e from https://github.com/rrwick/Badread.git)
#         run(f"badread simulate --seed {seed} --reference {Ref} --quantity 15x --length 15000,10000 --error_model pacbio --junk_reads 0 --random_reads 0 --chimeras 0 > {Reads} 2>{Data}/{id}_br.log.txt", gen_log, True, False)
    
#     aln_log = f"./{id}.log.txt"
#     long_gam = f"{Gams}/{id}_long.gam"
#     if not os.path.exists(long_gam):
#         run(f"{Bin} -x vg -f {Reads} -g {Graph} -a {long_gam}", aln_log)
    
#     for L, S, G in params:
#         clc_gam = f"{Gams}/{id}_clc_{L}_{S}_{G}.gam"
#         if not os.path.exists(clc_gam):
#             run(f"{Bin} --colinear-chaining -x vg -f {Reads} -g {Graph} -a {clc_gam} --colinear-gap {G} --colinear-split-len {L} --colinear-split-gap {S} ", aln_log)





# graph = sys.argv[1]
# reads = sys.argv[2]
# id = sys.argv[3]

# outdir = "out"
# if 4 < len(sys.argv):
# 	outdir = sys.argv[4]
# if not outdir.endswith(id + '/'):
# 	outdir = outdir + '/' + id + '/'

# def LoadGfaGraph(filename):
# 	VL, E = {}, {}
# 	for line in open(filename).readlines():
# 		if line[0] == 'S':
# 			# S 92533   A
# 			i, s = line[1:].strip().split()
# 			VL[int(i)] = s
# 		elif line[0] == 'L':
# 			# L 104890  +   104892  +   0M
# 			li, lr, ri, rr, ov = line[1:].strip().split()
# 			li, ri = int(li), int(ri)
# 			if li not in E:
# 				E[li] = []
# 			E[li].append(ri)
# 	return VL, E

# VL, E = LoadGfaGraph(graph)

# ed_global = lambda s1, s2 : edlib.align(s1, s2, mode='NW')['editDistance']
# ed_local = lambda s1, s2 : edlib.align(s1, s2, mode='HW')['editDistance']
# list2idx = lambda a: { a[i] : i for i in range(len(a)) }
# revc = lambda s: ''.join({"A":"T","T":"A","C":"G","G":"C"}[c] for c in s[::-1])

# import vg_pb2
# import gzip
# from google.protobuf.internal.encoder import _VarintBytes
# from google.protobuf.internal.decoder import _DecodeVarint32, _VarintDecoder
# def _VarintDecoder(mask):
# 	local_ord = ord
# 	def DecodeVarint(buffer, pos):
# 		result = 0
# 		shift = 0
# 		while 1:
# 			b = local_ord(buffer[pos])
# 			result |= ((b & 0x7f) << shift)
# 			pos += 1
# 			if not (b & 0x80):
# 				result &= mask
# 				return (result, pos)
# 			shift += 7
# 			# if shift >= 64:
# 			# 	raise _DecodeError('Too many bytes when decoding varint.')
# 	return DecodeVarint
# _DecodeVarint64 = _VarintDecoder((1 << 64) - 1)
# def read_alignments(gam_filename):
# 	with open(gam_filename, 'rb') as f:
# 		buf = gzip.GzipFile(fileobj=f).read()
# 		n = 0
# 		while n < len(buf):
# 			an, n = _DecodeVarint32(buf, n)
# 			for i in range(an):
# 				msg_len, n = _DecodeVarint32(buf, n)
# 				msg_buf = buf[n:n+msg_len]
# 				n += msg_len
# 				aln = vg_pb2.Alignment()
# 				aln.ParseFromString(msg_buf)
# 				yield aln

# def parse_alignment(aln):
# 	# bps = sum(len(VL[x.position.node_id]) for x in a.path.mapping)
# 	name = aln.name.split()[0]
# 	seq = ''

# 	rev_cnt = 0
# 	for x in aln.path.mapping:
# 		idx = x.position.node_id
# 		ll = VL[idx]
# 		if x.position.is_reverse:
# 			rev_cnt += 1
# 			seq += revc(ll)
# 		else:
# 			seq += ll
# 	return {'name':name, 'seq':seq, 'path_cnt':len(aln.path.mapping), 'revcnt':rev_cnt, 'path_bps':len(seq)}

# def parse_gam(filename):
# 	ret = {}
# 	for aln in read_alignments(filename):
# 		a = parse_alignment(aln)
# 		ret[a['name']] = a
# 	return ret

# seqs_long = parse_gam(f'{outdir}{id}_long.gam')
# seqs_clcs = parse_gam(f'{outdir}{id}_clc.gam')

# def read_fastq(fastq_filename):
# 	reads_lines = open(fastq_filename).readlines()
# 	for i, line in enumerate(reads_lines):
# 		if line[0] == '@':
# 			info = reads_lines[i].strip()
# 			# name = info.split()[0][1:]
# 			seq = reads_lines[i + 1].strip()
# 			yield (info, seq)

# seqs_read = {info.split()[0][1:] : (seq, info) for info, seq in read_fastq(f'{reads}')}

# class CSV:
# 	def __init__(self):
# 		self.h = []
# 		self.hidx = {}
# 		self.r = []
# 		self.ridx = {}
# 		self.data = []
# 	def add_headers(self, headers):
# 		for h in headers:
# 			if h not in self.hidx:
# 				self.h.append(h)
# 				self.hidx[h] = len(self.h) - 1
# 	def get_hids(self, headers):
# 		return [self.hidx[x] for x in headers]
# 	def add(self, row, hids = []):
# 		if len(hids) == 0:
# 			self.data.append(row[:])
# 		else:
# 			tmp = [''] * len(self.h)
# 			for i in range(len(hids)):
# 				tmp[hids[i]] = row[i]
# 			self.data.append(tmp[:])
# 	def save(self, filename):
# 		fout = open(filename, 'w')
# 		fout.write(','.join(self.h) + '\n')
# 		for d in self.data:
# 			if len(d) < len(self.h):
# 				d += [''] * (len(self.h) - len(d))
# 			fout.write(','.join(d) + '\n')
# 		fout.close()

# csv = CSV()
# csv.add_headers(['name', 'length', 'br_id_rate']) #0,1,2
# csv.add_headers(['long_pathcnt', 'long_path_bps', 'long_revcnt']) #3,4,5
# csv.add_headers(['clcs_pathcnt', 'clcs_path_bps', 'clcs_revcnt']) #6,7,8
# csv.add_headers(['long_align_rate']) #9
# csv.add_headers([
# 'global_ed_read_long',  #10
# # 'global_ed_long_true',
# 'global_ed_read_clcs',  #11
# # 'global_ed_clcs_true',
# # 'local_ed_read_long',
# # 'local_ed_long_read',
# # 'local_ed_true_long',
# # 'local_ed_long_true',
# # 'local_ed_read_clcs',
# # 'local_ed_clcs_read',
# # 'local_ed_true_clcs',
# # 'local_ed_clcs_true',
# # 'global_ed_read_clcs',
# # 'global_ed_long_clcs',
# # 'global_ed_read_true',
# # 'global_ed_clcs_true',
# ])

# reads_cnt = 0
# for name in seqs_read:
# 	reads_cnt += 1
# 	if reads_cnt % (len(seqs_read) // 5 + 1) == 0:
# 		print(reads_cnt, '/', len(seqs_read))
# 	seq, info = seqs_read[name]
# 	row = [''] * len(csv.h)
# 	row[0] = name
# 	for t in info.split():
# 		if t.startswith('length='):
# 			row[1] = t.split('=')[-1]
# 	row[2] = str('%.3f'%(float(info.split()[-1].split('=')[-1][:-1]) / 100))

# 	long_seq = ''
# 	if name in seqs_long:
# 		a = seqs_long[name]
# 		long_seq = a['seq']
# 		row[3] = str(a['path_cnt'])
# 		row[4] = str(a['path_bps'])
# 		row[5] = str(a['revcnt'])
# 		row[10] = str(ed_global(seq, long_seq))
# 	row[9] = str(len(long_seq) / len(seq))
# 	clcs_seq = ''
# 	if name in seqs_clcs:
# 		a = seqs_clcs[name]
# 		clcs_seq = a['seq']
# 		row[6] = str(a['path_cnt'])
# 		row[7] = str(a['path_bps'])
# 		row[8] = str(a['revcnt'])
# 		row[11] = str(ed_global(seq, clcs_seq))
# 	csv.add(row)

# csv.save(f'{outdir}{id}_summary.csv')
