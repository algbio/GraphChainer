#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <algorithm>
#include <random>

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

std::random_device rd;
std::mt19937 gen;
LL randint(LL l, LL r) {
	// random integer from [l .. r]
	return gen() % (r - l + 1) + l;
}

int main(int argc, char *argv[]) {
	// std::ios::sync_with_stdio(false);
	std::string input_filename = "", output_filename = "";
	LL anchors_count = 100;
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
					input_filename = peek(i + 1);
				else if (cmd == 'o')
					output_filename = peek(i + 1);
				else if (cmd == 'n')
					anchors_count = std::stoll(peek(i + 1));
			}
	}
	else {
		std::cout << "usage : ./a -f input.file -n number.of.anchors -o output.file";
	}
	loadGraph(input_filename, [&](LL i){}, [&](LL i, LL j){ add_edge(i, j); }, [&](LL N){ init(N); });
	std::cout << "Graph size " << N << " nodes, " << M << " edges" << std::endl;

	std::vector<Anchor> anchors(anchors_count);
	for (LL idx = 0; idx < anchors.size(); idx++) {
		Anchor &A = anchors[idx];
		LL length = randint(15, 30);
		for (LL i = randint(0, N - 1); A.path.size() < length; ) {
			A.path.push_back(i);
			LL d = 0;
			for (LL e = f[i]; e; e = t[e]) d++;
			if (d == 0) 
				break;
			d = randint(0, d - 1);
			LL e = f[i];
			while (d--) e = t[e];
			i = p[e];
		}
		A.x = randint(0, N - A.path.size());
		A.y = A.x + A.path.size() - 1;
	}

	std::cout << "generated " << anchors_count << " anchors, saving to " << output_filename << std::endl;

	saveAnchors(output_filename, anchors);
	return 0;
}