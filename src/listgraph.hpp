#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <random>

#include "Graph.hpp"
#include "Serializer.hpp"
#include "loader.hpp"

#include <lemon/lgf_reader.h>
#include <lemon/list_graph.h>

using namespace lemon;
using namespace std;

typedef long long LL;

// LL N, M;
std::vector<LL> f, p, t;
std::vector<std::string> lbs;
std::unordered_map<LL, LL> ids_map;

void init(LL N) {
	f.clear();
	f.resize(N, 0);
	t.clear();
	t.resize(1);
	p.clear();
	p.resize(1);
	lbs.clear();
	lbs.resize(N);
}

void add_edge(LL i, LL j) {
	p.push_back(j);
	t.push_back(f[i]);
	f[i] = t.size() - 1;
}

void loadGraphFile(const std::string &filename) {
	loadGraph(filename, [&](LL i, LL oid, std::string s){ lbs[i] = s; ids_map[oid] = i; }, [&](LL i, LL j){ add_edge(i, j); }, [&](LL N){ init(N); });
}