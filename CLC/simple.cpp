#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <algorithm>
// #include <iomanip>

#include "Graph.hpp"
#include "Serializer.hpp"
#include "loader.hpp"
#include "anchor.hpp"

#include <lemon/lgf_reader.h>
#include <lemon/list_graph.h>

using namespace lemon;
using namespace std;

typedef long long LL;

// LL N, M;
std::vector<LL> f, p, t;

void init(LL N) {
	f.clear();
	f.resize(N, 0);
	t.clear();
	t.resize(1);
	p.clear();
	p.resize(1);
}

void add_edge(LL i, LL j) {
	p.push_back(j);
	t.push_back(f[i]);
	f[i] = t.size() - 1;
}

bool reachable(LL S, LL T) {
	std::vector<LL> Q = { S }, d(N, -1);
	d[S] = 1;
	for (LL x = 0; x < Q.size() && d[T] == -1; x++) {
		LL u = Q[x];
		for (LL e = f[u]; e; e = t[e]) {
			LL v = p[e];
			if (d[v] == -1) {
				Q.push_back(v);
				d[v] = 1;
			}
		}
	}
	return d[T] != -1;
}

bool checkChain(const Chain &chain) {
	for (LL i = 0; i + 1 < chain.size(); i++) {
		LL S = chain[i].path.back(), T = chain[i + 1].path[0];
		if (chain[i].y >= chain[i + 1].y)
			return false;
		if (!reachable(S, T)) return false;
	}
	return true;
}
std::vector<LL> topo; // global topological order cache
std::vector<LL> toposort() {
	std::vector<LL> incd(N, 0);
	for (LL i = 0; i < N; i++)
		for (LL j = f[i]; j ; j = t[j])
			incd[p[j]]++;
	std::vector<LL> Q;
	for (LL i = 0; i < N; i++)
		if (incd[i] == 0)
			Q.push_back(i);
	for (LL i = 0; i < Q.size(); ) {
		LL s = Q[i++];
		for (LL j = f[s]; j ; j = t[j]) {
			LL t = p[j];
			incd[t]--;
			if (incd[t] == 0)
				Q.push_back(t);
		}
	}
	return Q;
}

typedef std::vector<LL> Path;
typedef std::vector<Path> PathCover;

PathCover greedyCover() {
	PathCover ret;
	std::vector<LL> covered(N, 0);
	LL covered_cnt = 0;
	std::vector<LL> incd(N, 0), pre(N, -1), dis(N, 0), Q(N, 0);
	while (covered_cnt < N) {
		for (LL i = 0; i < N; i++) {
			pre[i] = -1;
			dis[i] = 0;
			for (LL j = f[i]; j ; j = t[j])
				incd[p[j]]++;
		}
		LL Qsize = 0;
		for (LL i = 0; i < N; i++)
			if (incd[i] == 0)
				Q[Qsize++] = i;
		LL best = 0, best_i = -1;
		for (LL i = 0; i < Qsize; ) {
			LL s = Q[i++];
			if (covered[s] == 0)
				dis[s]++;
			if (best < dis[s]) {
				best = dis[s];
				best_i = s;
			}
			for (LL j = f[s]; j ; j = t[j]) {
				LL t = p[j];
				incd[t]--;
				if (dis[t] < dis[s]) {
					dis[t] = dis[s];
					pre[t] = s;
				}
				if (incd[t] == 0)
					Q[Qsize++] = t;
			}
		}
		Path tmp, path;
		for (LL i = best_i; i != -1; i = pre[i])
			tmp.push_back(i);
		std::reverse(tmp.begin(), tmp.end());
		LL l = 0, r = tmp.size() - 1;
		while (covered[tmp[l]]) l++;
		while (covered[tmp[r]]) r--;
		LL new_covered = 0;
		for (LL i = l; i <= r; i++) {
			path.push_back(tmp[i]);
			if (covered[tmp[i]] == 0)
				new_covered++;
			covered[tmp[i]]++;
		}
		covered_cnt += new_covered;
		// std::cout << "path #" << ret.size() << " : " << path.size() << " " << new_covered << " " << (N - covered_cnt) << std::endl;
		ret.push_back(path);
	}
	return ret;
}

struct flowGraph {
	LL N, M, S, T;
	std::vector<LL> f, p, t, c;

	flowGraph(LL NN) : N(NN+2) {
		init(N);
		S = NN;
		T = NN + 1;
	}

	void init(LL N) {
		f.clear();
		f.resize(N, 0);
		t.clear();
		t.resize(2);
		p = t;
		c = t;
	}

	void add_edge(LL i, LL j, LL cap) {
		// std::cout << "add edge " << i << " " << j << "  " << cap << std::endl;
		p.push_back(j);
		t.push_back(f[i]);
		c.push_back(cap);
		f[i] = t.size() - 1;
	}
};

PathCover shrink(PathCover &pc) {
	// graph should be DAG
	PathCover ret;
	LL K = pc.size(), inf = pc.size();
	std::vector<LL> covered(N, 0), starts(N, 0), ends(N, 0);
	std::map<std::pair<LL, LL>, LL> edge_covered;
	for (auto path : pc) {
		for (LL i = 0; i < path.size(); i++) {
			covered[path[i]]++;
			if (i > 0)
				edge_covered[{ path[i - 1], path[i] }]++;
		}
		starts[path[0]]++;
		ends[path.back()]++;
	}
	flowGraph fg(N * 2);
	// i_in = i, i_out = i + N
	// add r(i, j) = c(j,i) + f(i,j) - l(i,j)
	auto add = [&](LL i, LL j, LL cap, LL l, LL ff) {
		// std::cout << "add edge " << i << " " << j << "  " << cap << "  " << l << " " <<ff << std::endl;
		fg.add_edge(i, j, 0 + ff - l);
		fg.add_edge(j, i, cap - ff);
	};
	for (LL i = 0; i < N; i++)
		for (LL e = f[i]; e; e = t[e]) {
			// i_out -> j_in : [0, inf]
			LL j = p[e];
			LL ff = edge_covered.count({i, j}) ? edge_covered[{i, j}] : 0;
			add(i + N, j, inf, 0, ff);
		}
	for (LL i = 0; i < N; i++) {
		add(i, i + N, inf, 1, covered[i]);
		add(fg.S, i, inf, 0, starts[i]);
		add(i + N, fg.T, inf, 0, ends[i]);
	}
	
	LL total = inf;
	std::cout << "start shrinkink from width " << total << std::endl;
	std::vector<LL> Q(fg.N, 0), pre(fg.N, -1), d(fg.N, 0);
	while (1) {
		LL Qsize = 0;
		Q[Qsize++] = fg.S;
		for (LL i = 0; i < fg.N; i++) {
			pre[i] = -1;
			d[i] = 0;
		}
		d[fg.S] = 1;
		for (LL idx = 0; idx < Qsize && d[fg.T] == 0;) {
			LL i = Q[idx++];
			for (LL e = fg.f[i]; e; e = fg.t[e]) {
				LL j = fg.p[e]; 
				if (fg.c[e] > 0 && d[j] == 0) {
					d[j] = 1;
					pre[j] = e;
					Q[Qsize++] = j;
				}
			}
		}
		if (d[fg.T] == 0) break;
		Path tmp;
		LL flow = fg.c[pre[fg.T]];
		for (LL i = fg.T; ;) {
			tmp.push_back(i);
			LL e = pre[i];
			if (e == -1) break;
			flow = std::min(flow, fg.c[e]);
			i = fg.p[e ^ 1];
		}
		for (LL i = fg.T; ;) {
			LL e = pre[i];
			if (e == -1) break;
			fg.c[e] -= flow;
			fg.c[e^1] += flow;
			i = fg.p[e ^ 1];
		}
		if (flow == 0) exit(1);
		total -= flow;
		// std::cout << "Now shrink by " << flow << " to " << total << std::endl;
	}
	// convert flow back to path cover
	// ret = pc;
	// ret.resize(total);
	for (LL itr = 0; itr < total; itr++) {
		Path tmp;
		for (LL i = fg.S; i != fg.T; ) {
			if (0 <= i && i < N)
				tmp.push_back(i);
			LL nxt = -1;
			for (LL e = fg.f[i]; e; e = fg.t[e]) {
				LL j = fg.p[e];
				LL ff = fg.c[e] + ((i < N && i + N == j) ? 1 : 0);
				if ((e & 1) == 0 && ff > 0) {
					nxt = j;
					fg.c[e]--;
					break;
				}
			}
			if (nxt == -1) {
				std::cout << i << " not found nxt" << std::endl;
				return ret;
			}
			i = nxt;
		}
		ret.push_back(tmp);
	}
	return ret;
}

std::vector<std::vector<std::pair<LL, LL>>> forwards;
std::vector<std::vector<LL>> last2reach, paths;
void computeMPCIndex(const PathCover &pc) {
	LL K = pc.size();
	forwards.clear();
	forwards.resize(N);
	last2reach.clear();
	last2reach.resize(N, std::vector<LL>(K, -1));
	paths.clear();
	paths.resize(N);
	for (LL i = 0; i < K; i++)
		for (LL j = 0; j < pc[i].size(); j++) {
			last2reach[pc[i][j]][i] = j;
			paths[pc[i][j]].push_back(i);
		}
	for (LL i : topo) {
		for (LL e = f[i]; e; e = t[e]) {
			LL j = p[e];
			for (LL k = 0; k < K; k++)
				last2reach[j][k] = std::max(last2reach[j][k], last2reach[i][k]);
		}
	}
	for (LL i = 0; i < N; i++)
		for (LL k = 0; k < K; k++) {
			LL &idx = last2reach[i][k];
			if (idx != -1 && pc[k][idx] == i)
				idx--;
			if (idx != -1) {
				idx = pc[k][idx];
				forwards[idx].push_back({ i, k });
			}
		}
}

// struct TreeArray {
// 	LL n;
// 	vector<LL> a;
// 	TreeArray() {}
// 	TreeArray(LL n) {
// 		init(n);
// 	}
// 	void init(LL n) {
// 		this->n = n;
// 		a.clear();
// 		a.resize(n + 1, 0);
// 	}
// 	void add(LL i, LL x) {
// 		// a[0] += x;
// 		for (; i <= n; i += i & -i)
// 			a[i] = std::max(a[i], x);
// 	}
// 	// return sum{ a[1..i] }
// 	LL get(LL i) {
// 		LL ret = 0;
// 		for (; i > 0; i -= i & -i)
// 			ret = std::max(ret, a[i]);
// 		return ret;
// 	}
// 	// LL size() { return a[0]; }
// };

template<typename T>
struct SegmentTree {
	LL N;
	std::vector<T> t;
	LL ql, qr;
	T qx;

	LL get_upper_bit(LL N) {
		LL x = N & -N;
		while (x != N) {
			N ^= x;
			x = N & -N;
		}
		return x << 1;
	}
	static LL ID(LL l, LL r) {
		return l + r - 1;
	}
	SegmentTree(LL N_, T value) : N(get_upper_bit(N_)) {
		reset(value);
	}
	void reset(T value) {
		t.clear();
		t.resize(N*2, value);
	}
	void add(LL x, T value) {
		ql = x; qx = value;
		add_(0, N);
	}
	void add_(LL l, LL r) {
		t[ID(l, r)] = std::max(t[ID(l, r)], qx);
		if (l + 1 == r)
			return;
		LL m = (l + r) >> 1;
		if (ql < m)
			add_(l, m);
		else
			add_(m, r);
	}
	T RMQ(LL l, LL r) {
		if (l > r)
			return t[0];
		ql = l; qr = r + 1;
		return RMQ_(0, N);
	}
	T RMQ_(LL l, LL r) {
		// if (ql==27579&&qr==27605+1) cout<<"now "<<ID(l,r)<<"  "<<l<<" "<<r<<" : "<<t[ID(l,r)].first<<" "<<t[ID(l,r)].second<<endl;
		if (ql <= l && r <= qr)
			return t[ID(l, r)];
		LL m = (l + r) >> 1;
		if (qr <= m)
			return RMQ_(l, m);
		if (m <= ql)
			return RMQ_(m, r);
		return std::max(RMQ_(l, m), RMQ_(m, r));
	}
};

std::unordered_map<LL, LL> getSortedMap(std::vector<LL> a) {
	std::sort(a.begin(), a.end());
	a.erase(std::unique(a.begin(), a.end()), a.end());
	std::unordered_map<LL, LL> ret;
	for (LL i = 0; i < a.size(); i++)
		ret[a[i]] = i;
	return ret;
}

// input: pathcover, grpah, anchors
Chain sparseDP(const PathCover &pc, const std::vector<Anchor> &A) {
	// computeMPCIndex(pc);
	LL K = pc.size();
	std::pair<LL, LL> defaul_value = { -N*2, -1 };
	std::cout <<"defaul_value "<<defaul_value.first<<endl;
	std::vector<SegmentTree<std::pair<LL, LL>>> T(K, SegmentTree(N, defaul_value)), I(K, SegmentTree(N, defaul_value));
	std::vector<std::vector<LL>> starts(N), ends(N);
	std::vector<std::pair<LL, LL>> C(A.size());
	for (LL j = 0; j < A.size(); j++) {
		starts[A[j].path[0]].push_back(j);
		ends[A[j].path.back()].push_back(j);
		C[j] = { A[j].y - A[j].x + 1, -1 };
	}
	std::cout << "start dp" << std::endl;
	for (LL vidx = 0; vidx < N; vidx++) {
		LL v = topo[vidx];
		// std::cout << "now " << vidx << " : " <<v << " : " << starts[v].size() << "  " << ends[v].size() << std::endl;
		if (!starts[v].empty() && !ends[v].empty()) {
			std::vector<LL> ids = starts[v];
			for (LL j : ends[v])
				ids.push_back(j);
			std::sort(ids.begin(), ids.end(), [&](LL i, LL j) { 
				if (A[i].y != A[j].y)
					return A[i].y < A[j].y;
				return A[i].x < A[j].x;
			});
			ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
			// get local id count
			std::vector<LL> pos = { 0 };
			for (LL j : ids) {
				pos.push_back(A[j].x - 1);
				pos.push_back(A[j].x);
				pos.push_back(A[j].y);
			}
		// std::cout<<"step 1"<<endl;
			auto id_map = getSortedMap(pos);
			LL Size = (LL)id_map.size();
			SegmentTree tmpT(Size, defaul_value), tmpI(Size, defaul_value);
			for (LL j : ids) {
				if (A[j].path[0] == v) {
					// if (vidx ==13103)cout<<"n?? 0 " << "  "<<id_map[0]<<" "<< id_map[A[j].x - 1]<< endl;
					std::pair<LL, LL> q = tmpT.RMQ(id_map[0], id_map[A[j].x - 1]);
					// if (vidx ==13103)cout<<"n?? 1 " << endl;
				// if (j==44||j==9)cout<<v<<" _1_ "<<j<<" : " << A[j].x<<","<<A[j].y<<"  :"<<C[j].first<<" "<<C[j].second<<" <- "<<q.first<<" "<<q.second<<endl; 
					C[j] = std::max(C[j], {A[j].y - A[j].x + 1 + q.first, q.second});
					q = tmpI.RMQ(id_map[A[j].x], id_map[A[j].y]);
					// if (vidx ==13103)cout<<"n?? 2 " << endl;
				// if (j==44||j==9)cout<<v<<" _2_ "<<j<<" : " << A[j].x<<","<<A[j].y<<"  :"<<C[j].first<<" "<<C[j].second<<" <- "<<q.first<<" "<<q.second<<endl; 
					C[j] = std::max(C[j], {A[j].y + q.first, q.second});
					// if (vidx ==13103)cout<<"n?? 3 " << endl;
				}
				if (A[j].path.back() == v) {
					tmpT.add(id_map[A[j].y], {C[j].first, j});
					// if (vidx ==13103)cout<<"n?? 4 " << endl;
					// if (j==44||j==9)cout<<v<<" _2_ "<<j<<" : "<<"  add "<<C[j].first-A[j].y<<endl;
					tmpI.add(id_map[A[j].y], {C[j].first - A[j].y, j});
					// if (vidx ==13103)cout<<"n?? 5 " << endl;
				}
			}
		}
		// std::cout<<"step 2"<<endl;
		for (LL j : ends[v]) {
			for (LL k : paths[v]) {
				T[k].add(A[j].y, {C[j].first, j});
				I[k].add(A[j].y, {C[j].first - A[j].y, j});
					// if (j==44||j==9)cout<<v<<" add "<<j<<" : "<<"  add "<<A[j].y<<"  "<<C[j].first-A[j].y<<endl;
			}
		}
		// std::cout<<"step 3"<<endl;
		for (auto ui : forwards[v]) {
			LL u = ui.first, k = ui.second;
			for (LL j : starts[u]) {
				std::pair<LL, LL> q = T[k].RMQ(0, A[j].x - 1);
				// if (j==44||j==9)cout<<v<<" _3_ "<<j<<" : " << A[j].x<<","<<A[j].y<<"  :"<<C[j].first<<" "<<C[j].second<<" <- "<<q.first<<" "<<q.second<<endl;
				C[j] = std::max(C[j], {A[j].y - A[j].x + 1 + q.first, q.second}); 
				q = I[k].RMQ(A[j].x, A[j].y);
				// if (j==44||j==9)cout<<v<<" _4_ "<<j<<" : " << A[j].x<<","<<A[j].y<<"  :"<<C[j].first<<" "<<C[j].second<<" <- "<<q.first<<" "<<q.second<<endl;
				C[j] = std::max(C[j], {A[j].y + q.first, q.second});
			}
		}
		// std::cout<<"step 4"<<endl;
	}
	std::pair<LL, LL> best = {0, -1};
	for (LL j = 0; j < A.size(); j++) 
		best = std::max(best, { C[j].first, j });
	std::cout << "optimal coverage : " << best.first << std::endl;
	Chain ret;
	for (LL i = best.second; i != -1; i = C[i].second) {
		ret.push_back(A[i]);
		// std::cout << "now " << i << "  " << C[i].first << " "  << C[i].second << endl;
	}
	// for (LL i:std::vector<LL>({
	// 	96,15,75,62,54,68,30,59,90,61,21,27,57,12,78,83,29,66,7,99,53,52,3,24,25,47,5,33,85,49,74,88,10,93,94,65,50,31,36,56,22,43,19,82,95,97,70,69,9,44,84
	// 	})) {
	// 	std::cout << "now " << i << ": " << C[i].first << " " << C[i].second << "   " << A[i].path[0]<<"->"<<A[i].path.back()<<" ["<<A[i].x <<","<<A[i].y<<"]" << endl;
	// }
	std::reverse(ret.begin(), ret.end());
	return ret;
}

bool checkPathCover(const PathCover &pc) {
	std::vector<LL> covered(N, 0);
	for (auto path : pc) {
		for (LL i = 0; i < path.size(); i++) {
			covered[path[i]]++;
			if (i > 0) {
				bool found = false;
				for (LL e = f[path[i - 1]]; e && !found; e = t[e])
					found |= p[e] == path[i];
				if (!found) {
					std::cout << "edge not found : " << path[i - 1] << " -> " << path[i] << std::endl;
					return false;
				}
			}
		}
	}
	for (LL i = 0; i < N; i++)
		if (covered[i] == 0) {
			std::cout << "node not covered : " << i << std::endl;
			return false;
		}
	return true;
}



int main(int argc, char *argv[]) {
	// std::ios::sync_with_stdio(false);
	std::string graph_filename = "", anchor_filename = "";
	if (argc > 1) {
		auto peek = [&](int idx) {
			if (idx < argc)
				return std::string(argv[idx]);
			else {
				std::cerr << "error in augment " << idx << std::endl;
				exit(2);
			}
		};
		for (int i = 1; i < argc; i++)
			if (argv[i][0] == '-') {
				char cmd = argv[i][1];
				if (cmd == 'f')
					graph_filename = peek(i + 1);
				else if (cmd == 'a')
					anchor_filename = peek(i + 1);
			}
	}
	else {
		std::cout << "usage : ./a -f graph.input.file -a anchor.input.file";
	}
	loadGraph(graph_filename, [&](LL i){}, [&](LL i, LL j){ add_edge(i, j); }, [&](LL N){ init(N); });
	std::cout << "Graph size " << N << " nodes, " << M << " edges" << std::endl;
	
	std::vector<Anchor> anchors;
	if (anchor_filename != "") anchors = loadAnchors(anchor_filename);
	std::cout << "Loaded " << anchors.size() << " anchors" << std::endl;
	{
		LL sl = 0;
		std::vector<LL> ls;
		for (const Anchor &A : anchors) {
			ls.push_back(A.y - A.x + 1);
			sl += ls.back();
		}
		std::sort(ls.begin(), ls.end());
		double avg = sl / (double)anchors.size();
		std::cout << "avg anchor length " << avg << " from range ["<<ls[0]<<","<<ls.back()<<"] median " << ls[ls.size()/2] <<  endl;
	}

	// check DAG
	topo = toposort();
	std::cout << "topo size " << topo.size() << " -> " << (topo.size() == N ? "the graph is a DAG" : "the graph has cycles") << std::endl;
	if (topo.size() != N) {
		std::cerr << "not a DAG" << std::endl;
		return 1;
	}

	// greedy MPC
	PathCover gc = greedyCover();
	// {
	// 	auto s = Serializer("../chr22.greedycover.bin");
	// 	s << gc.size();
	// 	for (auto p : gc) s << p;
	// }
	// PathCover gc;
	// {
	// 	std::cout << "start loading greedy cover from cache file" << std::endl;
	// 	auto s = Deserializer("../chr22.greedycover.bin");
	// 	LL k;
	// 	s >> k;
	// 	for (LL i = 0; i < k; i++) {
	// 		Path p;
	// 		s >> p;
	// 		gc.push_back(p);
	// 	}
	// }
	// {
	// 	gc.push_back({ 0 });
	// 	gc.push_back({ 1 });
	// }
	std::cout << "greedy cover width " << gc.size() << std::endl;
	if (!checkPathCover(gc)) std::cout << "error : not a path cover!!" << std::endl;

	PathCover opt = shrink(gc);
	std::cout << "optimal width " << opt.size() << std::endl;
	if (!checkPathCover(opt)) std::cout << "error : not a path cover!!" << std::endl;

	computeMPCIndex(gc);
	std::cout << "MPCIndex is built" << std::endl;
	
	Chain chain = sparseDP(gc, anchors);
	if (!checkChain(chain))
		std::cout<<"error : not a chain" << std::endl;
	std::cout << "optimal chain coverage " << coverage(chain) << " chain size " << chain.size() << std::endl;


	return 0;
}