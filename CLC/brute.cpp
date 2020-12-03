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

Chain bruteCLC(const std::vector<Anchor> &anchors) {
	LL M = anchors.size();
	std::vector<LL> ids;
	std::vector<std::pair<LL, LL>> C(M);
	std::pair<LL, LL> best;
	std::vector<std::vector<LL>> starts(N);
	for (LL i = 0; i < M; i++) {
		ids.push_back(i);
		starts[anchors[i].path[0]].push_back(i);
	}
	std::sort(ids.begin(), ids.end(), [&](LL i, LL j){
		if (anchors[i].y != anchors[j].y)
			return anchors[i].y < anchors[j].y;
		return anchors[i].x < anchors[j].x;
	});
	std::vector<LL> Q(N), d(N, 0);
	LL Qsize = 0, flag = 0;
	for (LL i : ids) {
		const Anchor &A = anchors[i];
		C[i] = std::max(C[i], { A.y-A.x+1, -1 });
		Qsize = 0; flag++;
		Q[Qsize++] = A.path.back();
		d[Q[0]] = flag;
		for (LL idx = 0; idx < Qsize; idx++) {
			LL u = Q[idx];
			for (LL j : starts[u]) {
				if (anchors[j].y <= anchors[i].y)
					continue;
				LL value = C[i].first + anchors[j].y - std::max(anchors[j].x - 1, anchors[i].y);
				if (value > C[j].first)
					C[j] = { value , i };
			}
			for (LL e = f[u]; e; e = t[e]) {
				LL v = p[e];
				if (d[v] != flag) {
					d[v] = flag;
					Q[Qsize++] = v;
				}
			}
		}
		best = std::max(best, { C[i].first, i });
	}
	Chain ret;
	for (LL i = best.second; i != -1; i = C[i].second) {
		ret.push_back(anchors[i]);
		// std::cout << "now " << i << ": " << C[i].first << " " << C[i].second << "   " << anchors[i].path[0]<<"->"<<anchors[i].path.back()<<" ["<<anchors[i].x <<","<<anchors[i].y<<"]" << endl;
	}
	std::reverse(ret.begin(), ret.end());
	std::cout << "optimal coverage : " << best.first << std::endl;
	return ret;
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
	
	std::vector<Anchor> anchors = loadAnchors(anchor_filename);
	std::cout << "Loaded " << anchors.size() << " anchors" << std::endl;

	Chain chain = bruteCLC(anchors);
	if (!checkChain(chain))
		std::cout<<"error : not a chain" << std::endl;
	std::cout << "optimal chain coverage " << coverage(chain) << " chain size " << chain.size() << std::endl;

	return 0;
}