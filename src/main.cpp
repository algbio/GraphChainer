#include <iostream>
#include <fstream>
#include <sstream>
#include "Graph.hpp"
#include "Serializer.hpp"
#include "loader.hpp"

#include <lemon/lgf_reader.h>

#include <mpc/mpc.h>

using namespace lemon;

struct LemonGraph {
	ListDigraph g;
	// ListDigraph::NodeMap<int> label;
	LemonGraph(){}
};

void loadLemon(std::string filename) {
	ListDigraph g;
	digraphReader(g, filename).run();
	{
		int cnt = 0;
		for (ListDigraph::NodeIt n(g); n != INVALID; ++n)
			cnt++;
		std::cout << "Number of nodes: " << cnt << std::endl;
	}
	{
		int cnt = 0;
		for (ListDigraph::ArcIt a(g); a != INVALID; ++a)
			cnt++;
		std::cout << "Number of arcs: " << cnt << std::endl;
	}
}

void graphToLemon(Graph &G, ListDigraph &ret) {
	std::vector<ListDigraph::Node> xs;
	for (ID i = 0; i < G.VL.size(); i++) {
		xs.push_back(ret.addNode());
		// ret.label[xs.back()] = G.VL[i][0];
	}
	for (ID i = 0; i < G.VL.size(); i++)
		for (ID j : G.E[i]) {
			ret.addArc(xs[i], xs[j]);
		}
}


Graph loadTxtGraph(const std::string &filename) {
	// format specified in https://github.com/ParBLiSS/PaSGAL#graph-input-format as follow:
	// .txt is a simple human readable format. 
	// The first line indicates the count of total vertices (say n). 
	// Each subsequent line contains information of vertex i, 0 <= i < n. 
	// The information in a single line conveys its zero or more out-neighbor vertex ids, 
	// followed by its non-empty DNA sequence (either space or tab separated). 
	// For example, the following graph is a directed chain of four vertices: AC (id:0) -> GT (id:1) -> GCCGT (id:2) -> CT (id:3)
	Graph G;
	std::ifstream fin(filename);
	std::string s;
	LL N, n = -1;
	std::vector<std::pair<LL, LL>> edges;
	for (; std::getline(fin, s); n++) {
		if (n == -1)
			std::istringstream(s) >> N;
		else {
			LL i;
			char c = s.back();
			G.add_node(std::to_string(n), std::string(1, c));
			std::istringstream ssin(s.substr(0, s.length() - 1));
			while (ssin >> i)
				edges.push_back({ n, i });
		}
	}
	for (auto edge : edges)
		G.add_edge(std::to_string(edge.first), std::to_string(edge.second));
	return G;
}

Graph loadBinGraph(const std::string &filename) {
	Graph G;
	Deserializer(filename) >> G;
	return G;
}

Graph makeToyGraph(LL A = 100, LL B = 500) {
	std::cout << "building grid graph : " << A << " x " << B << std::endl;
	Graph G;
	for (LL i = 0; i <= A; i++)
		for (LL j = 0; j <= B; j++)
			G.add_node(std::to_string(i*(B+1)+j), "A");
	for (LL i = 0; i <= A; i++)
		for (LL j = 0; j <= B; j++) {
			if (i + 1 <= A)
				G.add_edge(std::to_string(i*(B+1)+j), std::to_string((i+1)*(B+1)+j));
			if (j + 1 <= B)
				G.add_edge(std::to_string(i*(B+1)+j), std::to_string(i*(B+1)+j+1));
		}
	return G;
}


int main(int argc, char *argv[]) {
	std::cout << "start loading" << std::endl;
	std::string filename = "";
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
					filename = peek(i + 1);
			}
	}
	else {
		std::cout << "usage : ./a -f input.file";
	}
	ListDigraph g;
	std::vector<ListDigraph::Node> xs;
	loadGraph(filename, 
		[&](LL i, std::string s){ xs.push_back(g.addNode()); }, 
		[&](LL i, LL j){ g.addArc(xs[i], xs[j]); }, 
		[&](LL N){ });

	// graphToLemon(G, g);

	std::cout << "MPC Min-Flow<MinCostFlow> reduction" << std::endl;
	auto mpc = MPC(g);
	std::cout << "path cover size : " << mpc.size() << std::endl;
	std::ofstream fout("mpc.txt");
	for (auto path: mpc) {
		for (auto vertex: path) {
			fout << g.id(vertex) << " , ";
		}
		fout << std::endl;
	}
	
	return 0;
}
