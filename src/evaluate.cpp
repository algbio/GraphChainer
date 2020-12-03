#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <unordered_set>
#include <unordered_map>
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
typedef std::vector<LL> Path;
typedef std::vector<Path> PathCover;

struct Reads {
	std::string name;
	std::string seq;
	std::string qs;
	static std::vector<Reads> LoadFastq(const std::string &filename) {
		std::vector<Reads> ret;
		std::ifstream fin(filename);
		std::string s;
		auto strip = [&](std::string &s) {
			while (!s.empty() && s.back() == '\n') s.pop_back();
		};
		while (std::getline(fin, s)) {
			Reads r;
			r.name = s.substr(1);
			strip(r.name);
			std::getline(fin, r.seq);
			strip(r.seq);
			std::getline(fin, s);
			std::getline(fin, r.qs);
			strip(r.qs);
			ret.push_back(r);
		}
		return ret;
	}
};

std::vector<std::string> splits_or(std::string s, std::string seps) {
	std::vector<std::string> ret = { "" };
	std::unordered_set<char> sep;
	for (char c : seps)
		sep.insert(c);
	for (char c : s)
		if (sep.count(c)) {
			if (!ret.back().empty())
				ret.push_back("");
		}
		else
			ret.back().push_back(c);
	if (ret.back().empty()) ret.pop_back();
	return ret;
}

bool reachable(LL S, LL T) {
	std::vector<LL> Q = { S }, d(N, -1);
	d[S] = 1;
	for (LL x = 0; x < Q.size() && d[T] == -1; x++) {
		LL u = Q[x];
		// std::cout << "? " << S << " "  << T << " : " << x << " : " << u << std::endl;
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

std::vector<LL> topo, topo_ids;
Path GetChainPath(const Chain &chain) {
	if (topo.empty()) {
		std::vector<LL> incd(N, 0);
		for (LL i = 0; i < N; i++)
			for (LL j = f[i]; j ; j = t[j])
				incd[p[j]]++;
		for (LL i = 0; i < N; i++)
			if (incd[i] == 0)
				topo.push_back(i);
		for (LL i = 0; i < topo.size(); ) {
			LL s = topo[i++];
			for (LL j = f[s]; j ; j = t[j]) {
				LL t = p[j];
				incd[t]--;
				if (incd[t] == 0)
					topo.push_back(t);
			}
		}
		topo_ids.resize(N, -1);
		for (LL i = 0; i < topo.size(); i++)
			topo_ids[topo[i]] = i;
	}
	std::unordered_set<LL> nodes;
	std::vector<std::pair<LL, LL>> ids;
	for (const Anchor &a : chain)
		for (LL i : a.path)
			if (!nodes.count(i)) {
				ids.push_back({ topo_ids[i], i });
				nodes.insert(i);
			}
	std::sort(ids.begin(), ids.end());
	ids.erase( std::unique( ids.begin(), ids.end() ), ids.end() );
	LL last = -1;
	Path ret;
	std::vector<LL> vis(N, 0), Q, pre(N);
	Q.reserve(N);
	LL flag = 1;
	for (auto pi : ids) {
		LL T = pi.second;
		if (last != -1) {
			Q.clear();
			flag++;
			Q.push_back(last);
			vis[last] = flag;
			for (LL i = 0; vis[T] != flag && i < Q.size(); ) {
				LL s = Q[i++];
				for (LL j = f[s]; j ; j = t[j]) {
					LL t = p[j];
					if (vis[t] != flag) {
						Q.push_back(t);
						vis[t] = flag;
						pre[t] = s;
					}
				}
			}
			std::vector<LL> tmp;
			for (LL i = T; i != last; i = pre[i])
				tmp.push_back(i);
			std::reverse(tmp.begin(), tmp.end());
			for (LL i = 0; i < tmp.size() - 1; i++)
				ret.push_back(tmp[i]);
		}
		ret.push_back(T);
		last = T;
	}
	return ret;
}

int main(int argc, char *argv[]) {
	// std::ios::sync_with_stdio(false);
	std::unordered_map<char, std::string> augs;
	augs['p'] = "path.bin";
	augs['o'] = "out.csv";
	if (argc > 1) {
		auto peek = [&](int idx) {
			if (idx < argc)
				return std::string(argv[idx]);
			else {
				std::cerr << "error no enough augment at position " << idx << std::endl;
				exit(2);
			}
		};
		for (int i = 1; i < argc; i++)
			if (argv[i][0] == '-')
				augs[argv[i][1]] = peek(i + 1);
	}
	if (argc <= 1 || augs['f'].empty()) {
		std::cout << "usage : ./evaluate -f graph.file -o output.csv.file -p path.bin -q out_chainpaths.fasta -c chains.txt -a anchors.txt -r reads.fastq";
		return 0;
	}
	loadGraphFile(augs['f']);

	Path p;
	Deserializer(augs['p']) >> p;
	std::vector<LL> ref_node_ids;
	for (LL i : p)
		for (char c : lbs[i]) 
			ref_node_ids.push_back(i);

	std::vector<Reads> reads = Reads::LoadFastq(augs['r']);
	std::unordered_map<std::string, LL> reads_ids_map;
	for (LL i = 0; i < reads.size(); i++) {
		std::string name = reads[i].name;
		//714d6421-ec9b-743e-964b-35b2abc85136 tmp,+strand,682500-686370 length=3743 error-free_length=3872 read_identity=82.13%
		name = name.substr(0, name.find(' '));
		reads_ids_map[name] = i;
	}

	std::vector<AnchorSet> anchor_sets = loadAnchorSets(augs['a']);
	std::unordered_map<std::string, LL> anchor_ids_map;
	for (LL i = 0; i < anchor_sets.size(); i++) {
		//714d6421-ec9b-743e-964b-35b2abc85136
		anchor_ids_map[anchor_sets[i].name] = i;
	}
	{ // map ids
		for (AnchorSet &as : anchor_sets) {
			for (Anchor &a : as.anchors)
				for (LL &i : a.path) {
					if (ids_map.count(i))
						i = ids_map[i];
					else {
						std::cout << "not found " << i << std::endl;
					}
				}
		}
	}

	std::unordered_map<std::string, std::vector<LL>> chains_map;
	{
		std::ifstream fin(augs['c']);
		std::string name, s;
		while (std::getline(fin, name)) {
			while (!name.empty() && name.back() == '\n') name.pop_back();
			std::getline(fin, s);
			std::istringstream ssin(s);
			std::vector<LL> &vec = chains_map[name];
			LL x;
			while (ssin >> x)
				vec.push_back(x);
		}
	}

	std::ofstream fout(augs['o']);
	std::ofstream fasta_out(augs['q']);
	fout << "id,length,error-free_length,start,end,path_length,short_reads,anchors_count,chain_nodes,chain_coverage,is_chain,overlap_nodes,overlap_raw_bps,chainpath_nodes,chainpath_bps,chainpath_overlap_nodes,chainpath_overlap_bps"  << std::endl;
	for (auto pa : reads_ids_map) {
		LL read_idx = pa.second;
		std::vector<std::string> tokens = splits_or(reads[read_idx].name, " ,=");
		std::string name = tokens[0];
		// 0                                    1   2       3             4      5    6                 7    8             9
		//714d6421-ec9b-743e-964b-35b2abc85136 tmp,+strand,682500-686370 length=3743 error-free_length=3872 read_identity=82.13%
		std::vector<std::string> lrs = splits_or(tokens[3], "-");
		LL l = std::stoll(lrs[0]), r = std::stoll(lrs[1]);
		LL len = std::stoll(tokens[5]);
		LL len_ef = std::stoll(tokens[7]);
		LL short_count = (len - 150 + 1) / 36 + 1; 
		LL path_length = 0, anchors_count = 0, chain_nodes = 0, chain_coverage = 0, overlap_nodes = 0, overlap_raw_bps = 0;
		std::string is_chain = "";
		LL chainpath_nodes = 0,chainpath_bps = 0,chainpath_overlap_nodes = 0, chainpath_overlap_bps = 0;
		std::string chainpath_seq = "";
		if (anchor_ids_map.count(name)) {
			Path tp; //true_path
			for (LL i = l; i <= r; i++)
				if (tp.empty() || ref_node_ids[i] != tp.back())
					tp.push_back(ref_node_ids[i]);
			path_length = tp.size();
			std::unordered_set<LL> tp_set(tp.begin(), tp.end());
			const AnchorSet &anchors = anchor_sets[anchor_ids_map[name]];
			anchors_count = anchors.anchors.size();
			if (chains_map.count(name)) {
				// std::cout<<"now " << name << " : " << anchor_ids_map.count(name) << " " << chains_map.count(name) <<std::endl;
				std::vector<LL> chain_ids = chains_map[name];
				Chain chain;
				for (LL i : chain_ids) {
					// std::cout << i << " " << anchors.anchors.size() << std::endl;
					chain.push_back(anchors.anchors[i]);
				}
		
				chain_coverage = coverage(chain);
				is_chain = checkChain(chain) ? "Yes" : "No";
				std::unordered_set<LL> overlaps, chain_nodes_set;
				for (const Anchor &a : chain)
					for (LL i : a.path) {
						if (tp_set.count(i)) 
							overlaps.insert(i);
						chain_nodes_set.insert(i);
					}
				chain_nodes = chain_nodes_set.size();
				overlap_nodes = overlaps.size();
				for (LL i : overlaps)
					overlap_raw_bps += lbs[i].length();
				Path chainpath = GetChainPath(chain);
				chainpath_nodes = chainpath.size();
				for (LL i : chainpath) {
					chainpath_bps += lbs[i].length();
					chainpath_seq += lbs[i];
				}
				overlaps.clear();
				for (LL i : chainpath)
					if (tp_set.count(i)) 
						overlaps.insert(i);
				chainpath_overlap_nodes = overlaps.size();
				for (LL i : overlaps)
					chainpath_overlap_bps += lbs[i].length();
			}
		}
		fout << name << "," << len << "," << len_ef << "," << l << "," << r << "," << path_length << "," << short_count << "," << anchors_count << "," << chain_nodes << "," << chain_coverage << "," << is_chain << "," << overlap_nodes << "," << overlap_raw_bps << "," << chainpath_nodes << "," << chainpath_bps << "," << chainpath_overlap_nodes << "," << chainpath_overlap_bps << std::endl;
		fasta_out << ">" << name << std::endl;
		fasta_out << chainpath_seq << std::endl;
	}
	fout.close();
	fasta_out.close();
	return 0;
}