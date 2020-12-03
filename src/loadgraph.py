
def loadgfa(filename):
    V, E, L = {}, set(), {}
    N, M = 0, 0
    for line in open(filename).readlines():
        if line[0] == 'S':
            c, id, label = line.split()
            V[int(id)] = N
            L[N] = label
            N += 1
        elif line[0] == 'L':
            c, lid, l, rid, r, o = line.split()
            lid, rid = V[int(lid)], V[int(rid)]
            E.add((lid, rid))
            M += 1
    return V, E, L


G = loadgfa('../data/LRC/LRC.gfa')
V, E, L = G
