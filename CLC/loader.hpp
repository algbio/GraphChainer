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

#include <lemon/lgf_reader.h>
#include <lemon/list_graph.h>

using namespace lemon;
using namespace std;

typedef long long LL;

LL N, M;

std::string getExtension(const std::string &s, const std::string &sep = ".") {
	size_t i = s.rfind(sep);
	if (i != std::string::npos)
		return s.substr(i + 1);
	return "";
}

template<typename FV, typename FE, typename FI>
void loadTxtGraph(const std::string &filename, FV add_node, FE add_edge, FI init) {
	// format specified in https://github.com/ParBLiSS/PaSGAL#graph-input-format as follow:
	// .txt is a simple human readable format. 
	// The first line indicates the count of total vertices (say n). 
	// Each subsequent line contains information of vertex i, 0 <= i < n. 
	// The information in a single line conveys its zero or more out-neighbor vertex ids, 
	// followed by its non-empty DNA sequence (either space or tab separated). 
	// For example, the following graph is a directed chain of four vertices: AC (id:0) -> GT (id:1) -> GCCGT (id:2) -> CT (id:3)
	std::ifstream fin(filename);
	std::string s;
	std::getline(fin, s);
	N = std::stoll(s);
	init(N);
	for (LL i = 0; i < N; i++)
		add_node(i);
	M = 0;
	std::vector<std::pair<LL, LL>> edges;
	for (LL n = 0; std::getline(fin, s); n++) {
		std::istringstream ssin(s.substr(0, s.length() - 1));
		for (LL i; ssin >> i; M++)
			edges.push_back({ n, i });
	}
	for (auto edge : edges)
		add_edge(edge.first, edge.second);
}


template<typename FV, typename FE, typename FI>
void loadBinGraph(const std::string &filename, FV add_node, FE add_edge, FI init) {
	Graph G;
	Deserializer(filename) >> G;
	N = G.VL.size();
	init(N);
	for (LL i = 0; i < N; i++)
		add_node(i);
	M = 0;
	for (LL i = 0; i < N; i++)
		for (LL j : G.E[i]) {
			add_edge(i, j);
			M++;
		}
}

template<typename FV, typename FE, typename FI>
void loadLemonGraph(const std::string &filename, FV add_node, FE add_edge, FI init) {
	ListDigraph g;
	digraphReader(g, filename).run();
	N = 0; M = 0;
	for (ListDigraph::NodeIt n(g); n != INVALID; ++n) {
		add_node(N);
		N++;
	}
	init(N);
	for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
		M++;
		add_edge(g.id(g.source(a)), g.id(g.target(a)));
	}	
}

template<typename FV, typename FE, typename FI>
void loadVgGraph(const std::string &filename, FV add_node, FE add_edge, FI init) {
}


template<typename FV, typename FE, typename FI>
void loadGfaGraph(const std::string &filename, FV add_node, FE add_edge, FI init) {
	std::string s;
	char c;
	std::unordered_map<LL, LL> ids;
	{ // count nodes
		std::ifstream fin(filename);
		while (std::getline(fin, s)) {
			istringstream ssin(s);
			ssin >> c;
			if (c == 'S') {
				string id, label;
				ssin >> id >> label;
				LL i = std::stoll(id);
				if (!ids.count(i))
					ids[i] = ids.size();
				add_node(ids[i]);
			}
		}
	}
	N = (LL)ids.size();
	init(N);
	M = 0;
	{ // add edges
		std::ifstream fin(filename);
		while (std::getline(fin, s)) {
			istringstream ssin(s);
			ssin >> c;
			if (c == 'L') {
				char l, r;
				string lid, rid, overlap;
				ssin >> lid >> l >> rid >> r >> overlap;
				LL li = ids[std::stoll(lid)], ri = ids[std::stoll(rid)];
				// if (overlap != "0M" || l != '+' || r != '+') {
				// 	cerr << "line " << n << " non-standard edge : " << s << endl;
				// 	return;
				// }
				add_edge(li, ri);
				M++;
			}
		}
	}
}

template<typename FV, typename FE, typename FI>
void loadGraph(const std::string &filename, FV fv, FE fe, FI init) {
	if (filename != "") {
		std::cout << "start loading from " << filename << std::endl;
		std::string ext = getExtension(filename);
		if (ext == "txt")
			loadTxtGraph(filename, fv, fe, init);
		else if (ext == "bin")
			loadBinGraph(filename, fv, fe, init);
		else if (ext == "vg")
			loadVgGraph(filename, fv, fe, init);
		else if (ext == "lgf")
			loadLemonGraph(filename, fv, fe, init);
		else if (ext == "gfa") 
			loadGfaGraph(filename, fv, fe, init);
		else {
			std::cout << "unknown input format: " << ext << ", trying loading as vg format" << std::endl;
			loadVgGraph(filename, fv, fe, init);
		}
	}
	else {

	}
}


// void makeToyGraph(LL A = 100, LL B = 500) {
// 	std::cout << "building grid graph : " << A << " x " << B << std::endl;
// 	N = (A+1)*(B+1);
// 	init((A+1)*(B+1));
// 	for (LL i = 0; i <= A; i++)
// 		for (LL j = 0; j <= B; j++) {
// 			if (i + 1 <= A)
// 				add_edge(i*(B+1)+j, (i+1)*(B+1)+j);
// 			if (j + 1 <= B)
// 				add_edge(i*(B+1)+j, i*(B+1)+j+1);
// 		}
// 	M = p.size() - 2;
// }
