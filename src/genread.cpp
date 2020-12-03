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
#include "listgraph.hpp"

#include <lemon/lgf_reader.h>
#include <lemon/list_graph.h>

using namespace lemon;
using namespace std;

typedef long long LL;

std::random_device rd;
std::mt19937 gen;
LL randint(LL l, LL r) {
	// random integer from [l .. r]
	return gen() % (r - l + 1) + l;
}

std::vector<LL> F;
LL getf(LL i) {
	if (F[i] != i)
		return F[i] = getf(F[i]);
	else
		return i;
}

typedef std::vector<LL> Path;
typedef std::vector<Path> PathCover;

Path randomPath(LL start) {
	Path ret;
	for (LL now = start; ; ) {
		ret.push_back(now);
		std::vector<LL> ids;
		for (LL e = f[now]; e; e = t[e])
			ids.push_back(p[e]);
		if (ids.empty())
			break;
		now = ids[randint(0, ids.size() - 1)];
	}
	return ret;
}

Path randomPath() {
	std::vector<LL> ID(N, 0), OD(N, 0);
	for (LL i = 0; i < N; i++)
		for (LL e = f[i]; e; e = t[e]) {
			LL j = p[e];
			ID[j]++; OD[i]++;
		}
	std::vector<LL> ss;
	for (LL i = 0; i < N; i++)
		if (ID[i] == 0)
			ss.push_back(i);
	return randomPath(ss[randint(0, ss.size() - 1)]);
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
			}
	}
	else {
		std::cout << "usage : ./a -f input.file -o output.file";
	}
	loadGraphFile(input_filename);
	std::cout << "Graph size " << N << " nodes, " << M << " edges" << std::endl;

	std::ofstream fout(output_filename);
	fout << ">tmp" << endl;
	Path p = randomPath();
	LL l = 0, LLEN = 100;
	// for (LL i= 0 ; i < p.size()&& i< 10; i++)cout << "path " << i << " : " << p[i] << " : " << lbs[i] << endl;
	for (LL i : p) {
		for (char c : lbs[i]) {
			fout << c;
			l++;
			if (l % LLEN == 0) fout << endl;
		}
		// LL ll = l;
		// l += lbs[i].length();
		// fout << lbs[i];
		// if (l / LLEN != ll / LLEN)
		// 	fout << endl;
	}
	fout << endl;
	fout.close();
	Serializer("path.bin") << p;
	std::cout << "path length " << p.size() << " nodes, " << l << " bp" << std::endl;
	return 0;
}