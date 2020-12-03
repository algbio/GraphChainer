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

outf = open('out.csv', 'w')
outf.write('Graph,#nodes,#edges,DAG?,lemon mpc,lemon mpc time,greedywidth,greedy time,width,total time,ret\n')

os.system('mkdir -p logs')
for filename in filenames:
    filen = filename.split('/')[-1]
    name = filen.split('.')[0]
    print('testing', name)
    t = time.perf_counter()
    os.system('timeout 40m ./lemon_gc -f %s > logs/log_lemon_%s.txt' % (filename, name))
    t1 = time.perf_counter() - t
    os.system('timeout 40m ./simple_gc -f %s > logs/log_simple_%s.txt' % (filename, name))
    t2 = time.perf_counter() - t - t1
    lemon_time = t1 / 60
    total_time = t2 / 60
    logl = open('logs/log_lemon_%s.txt' % (name)).read()
    lw = get_nxt_token(logl, 'path cover size : ')
    logs = open('logs/log_simple_%s.txt' % (name)).read()
    N = get_nxt_token(logs, 'Graph size ')
    M = get_nxt_token(logs, ' nodes, ')
    dag = "no" if "not a DAG" in logs else "yes"
    gw = get_nxt_token(logs, 'greedy cover width ')
    gwt = ''
    width = get_nxt_token(logs, 'optimal width ')
    ret = "Correct" if width == lw else "Wrong"
    outf.write(f'{name},{N},{M},{dag},{lw},{lemon_time},{gw},{gwt},{width},{total_time},{ret}\n')
    outf.flush()
    # break
