import os, time

# os.system('make gc')
# os.system('make main')

filenames = []
filenames.append('../data/LRC/LRC.gfa')
filenames.append('../data/MHC/MHC1.gfa')
for i in range(1, 21):
    filenames.append('../data/ariel_large_20/graph_%d.lgf' % i)
    
filenames.append('../data/MHC/MHC2.txt')
filenames.append('../data/chr22/chr22.graph.bin')

def get_nxt_token(s, ss):
    i = s.find(ss)
    if i == -1:
        return ''
    return s[i+len(ss):].split()[0]

outf = open('out_tests.csv', 'w')
outf.write('Graph,#nodes,#edges,width,#anchors,avg_anchor_length,l,r,median,brute_coverage,brute_size,brute_time,CLC_coverage,CLC_size,CLC_time,comment\n')

os.system('mkdir -p logs')
for filename in filenames:
    filen = filename.split('/')[-1]
    name = filen.split('.')[0]
    for N in [100,1000,10000]:
        print('testing', name, 'size', N)
        afn = filename[:-len(filen)] + name + ('_%d'%N) + '.anchors'
        os.system('./gen -f %s -n %d -o %s > logs/log_gen_%s_%d.txt' % (filename, N, afn, name, N));
        t = time.perf_counter()
        os.system('timeout 40m ./brute -f %s -a %s > logs/log_brute_%s_%d.txt' % (filename, afn, name, N))
        t1 = time.perf_counter() - t
        os.system('timeout 40m ./simple_gc -f %s -a %s > logs/log_simple_%s_%d.txt' % (filename, afn, name, N))
        t2 = time.perf_counter() - t - t1

        logb = open('logs/log_brute_%s_%d.txt' % (name, N)).read()
        logs = open('logs/log_simple_%s_%d.txt' % (name, N)).read()

        GN = get_nxt_token(logs, 'Graph size ')
        GM = get_nxt_token(logs, ' nodes, ')
        width = get_nxt_token(logs, 'optimal width ')
        aal = get_nxt_token(logs, 'avg anchor length ')
        alx = get_nxt_token(logs, ' from range ')
        aall, aalr = alx[1:-1].split(',')
        aalm = get_nxt_token(logs, 'median ')

        brute_c = get_nxt_token(logb, 'optimal chain coverage ')
        brute_s = get_nxt_token(logb, 'chain size ')
        brute_time = t1 / 60

        clc_c = get_nxt_token(logs, 'optimal chain coverage ')
        clc_s = get_nxt_token(logs, 'chain size ')
        clc_time = t2 / 60
        ret = "Correct" if brute_c==clc_c else "Wrong"
        if "error" in logs:
            ret = "Error"
        print(ret)
        if ret != "Correct":
            exit(1)

        outf.write(f'{name},{GN},{GM},{width},{N},{aal},{aall},{aalr},{aalm},{brute_c},{brute_s},{brute_time},{clc_c},{clc_s},{clc_time},{ret}\n')
        outf.flush()
