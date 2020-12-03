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

	{
		std::vector<LL> ID(N, 0), OD(N, 0);
		for (LL i = 0; i < N; i++)
			for (LL e = f[i]; e; e = t[e]) {
				LL j = p[e];
				ID[j]++; OD[i]++;
			}
		std::map<std::pair<LL, LL>, LL> ids;
		for (LL i = 0; i < N; i++)
			ids[{ ID[i], OD[i] }]++;
		LL o = 0;
		for (auto i : ids) {
			LL ii = i.first.first;
			LL io = i.first.second;
			if (ii == 0 || io == 0 || ii*io==1)
				std::cout << i.first.first << " " << i.first.second << " : " << i.second << std::endl;
			else o+=i.second;
		}
		cout << "others : " << o << endl;
	}
	
	{
		F.clear(); F.resize(N, -1);
		for (LL i = 0; i < N; i++)
			F[i] = i;
		for (LL i = 0; i < N; i++)
			for (LL e = f[i]; e; e = t[e]) {
				LL j = p[e];
				LL fi = getf(i), fj = getf(j);
					F[getf(i)] = getf(j);
			}
		LL count = 0;
		for (LL i = 0; i < N; i++)
			if (F[i] == i)
				count++;
		std::cout << "components count " << count << std::endl;
	}
	
	return 0;
}