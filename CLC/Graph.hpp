#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <algorithm>
#include <vector>
#include <string>
#include <queue>
#include <iostream>
#include <functional>
#include <unordered_set>
#include <unordered_map>


using namespace std;
typedef long long LL; // default int type
typedef long long ID; // default id type
// typedef string LABEL; // default label type

struct Graph {
public:
    /// Probably these can be converted to vector of integers if we map the node's ids to [0...|V|-1]
    unordered_map<string, ID> id_map; // map raw node id to int id
    vector<string> VL; // <id, label>
    vector<vector<ID>> E; // assume no overlap
    vector<LL> D; // in-degree

    struct Path {
        string name;
        vector<int> ids, overlaps;
    };
    vector<Path> Paths; // paths, such as the reference sequence
public:
    Graph() {}

	template<typename Ar>
	Ar& Serialize(Ar &ar) {
		return ar & id_map & VL & E;
	}

    void reset_size(LL size = -1) {
        size = max(size, (LL)id_map.size());
        if (size > min(VL.size(), E.size())) {
            VL.resize(size);
            E.resize(size);
            D.resize(size);
        }
    }

    // It creates a Graph from a *.gfa file
    Graph(const string &filename) {
        cout << "Loading from " << filename << endl;
        if (filename.substr(filename.size() - 4) != ".gfa") {
            cerr << "can only load .gfa files for now" << endl;
            return;
        }
        ifstream fin(filename);
        string s;
        for (LL n = 0; getline(fin, s); n++) {
            istringstream ssin(s);
            char c;
            ssin >> c;
            if (c == 'H')
                continue;
            else if (c == 'S') {
                string id, label;
                ssin >> id >> label;
                if (!add_node(id, label)) {
                    cerr << "line " << n << " error when adding node : " << s << endl;
                    return;
                }
            }
            else if (c == 'P') {
                Path path;
                string ids, overlaps, token;
                ssin >> path.name >> ids >> overlaps;
                istringstream idsin(ids), overlapsin(overlaps);
                while (getline(idsin, token, ',')) {
                    // assert(token[token.size() - 1] == '+');
                    if (token[token.size() - 1] != '+') {
                        cerr << "unsupported path with '-' sign" << endl;
                        break;
                    }
                    token = token.substr(0, token.size() - 1);
                    path.ids.push_back(get_id(token));
                }
                while (getline(overlapsin, token, ',')) {
                    // assert(token[token.size() - 1] == 'M');
                    path.overlaps.push_back(stoll(token.substr(0, token.size() - 1)));
                }
                Paths.push_back(path);
                cout << "path " << Paths.size() << " : " << path.ids.size() << " nodes" << endl;
            }
            else if (c == 'L') {
                char l, r;
                string lid, rid, overlap;
                ssin >> lid >> l >> rid >> r >> overlap;
                if (overlap != "0M" || l != '+' || r != '+' || !add_edge(lid, rid)) {
                    cerr << "line " << n << " non-standard edge : " << s << endl;
                    return;
                }
            }
            else {
                cerr << "line " << n << " unknown format : " << s.substr(0, 100) << endl;
                return;
            }
        }
    }

    ID get_id(const string &id) {
        ID i;
        if (!id_map.count(id)) {
            i = ID(id_map.size());
            id_map[id] = i;
        }
        else 
            i = id_map[id];
        return i;
    };

    // Add a node with its respective label only if is a new node_id, or old node_id with same label
    bool add_node(const string &id, const string &label) {
        ID i = get_id(id);
        if (i == VL.size()) {
            VL.push_back(label);
            return true;
        }
        else {
            if (i > VL.size())
                cerr << "node idx out of range" << endl;
            else if (i < VL.size() && VL[i] != label)
                cerr << "ignored duplicate node id with different labels! " << id << " : " << VL[i] << " <- " << label << endl;
            return false;
        }
    }

    // Add edge from node with label l to node with label r
    // and update in-degree of r
    /// What if I call this when the edge is alread there? -> duplicated in E, double counted in D
    bool add_edge(const string &l, const string &r) {
        ID li = get_id(l), ri = get_id(r); 
        reset_size();
        E[li].push_back(ri);
        D[ri]++;
        return true;
    }

    // It recomputes the in-degree based on E
    void reset_degree() {
        D.clear();
        D.resize(VL.size(), 0);
        for (const vector<ID> &vec : E)
            for (const ID &r : vec)
                D[r]++;
    }

    // It prints the number of vertices,
    // the length of the concatenation of the labels
    // and the number of edges of the graph
    void print_statistics() const {
        cout << VL.size() << " vertices" << endl;

        LL total = 0;
        unordered_set<char> alphabet;
        for (const auto &p : VL) {
            total += p.length();
            unordered_set<char> tmp(p.begin(), p.end());
            alphabet.insert(tmp.begin(), tmp.end());
        }
        cout << total << " bases" << endl;
        string alphabet_s(alphabet.begin(), alphabet.end());
        sort(alphabet_s.begin(), alphabet_s.end());
        cout << "alphabet size " << alphabet.size() << " : " << alphabet_s << endl;

        LL edgecnt = 0;
        for (const vector<ID> &vec : E)
            edgecnt += vec.size();
        cout << edgecnt << " edges" << endl;
    }

    LL get_components() const {
        vector<LL> F(VL.size());
        for (ID i = 0; i < VL.size(); i++) 
            F[i] = i;
        std::function<ID(ID)> getf;
        getf = [&F, &getf] (ID x) -> ID {
            return (F[F[x]] == F[x]) ? F[x] : (F[x] = getf(F[x]));
        };
        for (ID s = 0; s < VL.size(); s++)
            for (ID t : E[s]) {
            	ID fs = getf(s), ft = getf(t);
                if (getf(s) != getf(t))
                	if (rand()&1)
	                    F[F[s]] = F[t];
	               	else
	               		F[F[t]] = F[s];
	        }
        LL count = 0;
        for (ID s = 0; s < VL.size(); s++)
            if (F[s] == s)
                count++;
        return count;
    }


    template<typename TV, typename TE>
    void topo_sort_fn(TV fv, TE fe) {
        reset_degree();
        queue<ID> Q;
        // If a vertex is a source (no incoming edges) we apply fv and push it to the queue
        /// Probably we will only have a couple of sources, so maybe no need
        /// to check all the vertices
        for (ID i = 0; i < VL.size(); i++) 
            if (D[i] == 0) {
                fv(i);
                Q.push(i);
            }
        // Queue-based algorithm for topological-sort
        // first process the sources, then mark the edges you visit
        // if a node has all its incoming edges marked, add it to the queue
        while (!Q.empty()) {
            ID s = Q.front();
            Q.pop();
            for (ID t : E[s]) {
                D[t]--;
                fe(s, t);
                if (D[t] == 0) {
                    fv(t);
                    Q.push(t);
                }
            }
        }
    }

    // // Use topological sort to count the number of vertices
    // // in a topological-sort traversal of the graph
    // // if this number is equal to the number of vertices
    // // then the graph is a DAG
    LL topo_sort_visit() {
        LL visit = 0;
        topo_sort_fn([&](const ID &s){ visit++; }, [&](const ID &s, const ID &t){ });
        return visit;
    }

    // // It computes and outputs info about a path cover of the extended graph,
    // // That's it the graph where a node with a label is extended as a path
    // // with the character labels.
    // // This not necessarily minimum, but at least a O(log(|V|) approximation
    // // It is computed by taking the path with most uncovered nodes each time.
    // // To do this computation, the graph is traversed in topological order.
    // // When processing a sink we take the max size ending at its in-neighbours
    // // and add the label length of the nodes if it is still uncovered.
    // //
    // // At the end of each path computation we output:
    // // - Number of characters covered by the path
    // // - Number of uncovered nodes (string labeled) in the path
    // // - Number of nodes (string labeled) in the path
    // // - Id of the last node in the path
    // // - Total number of nodes (string labeled) by this and all previous paths)
    // // - Number of nodes (string labeled) in the graph
    // // only first $limit paths; will early stop with incomplete cover if limit reached
    vector<vector<ID>> greedy_MPC(LL limit = -1) {
        vector<LL> covered(VL.size(), 0), dis;
        vector<ID> pre;
        LL cover_cnt = 0;
        vector<vector<ID>> paths;
        while (cover_cnt < VL.size()) {
            dis.clear();
            dis.resize(VL.size(), 0);
            pre.clear();
            pre.resize(VL.size(), -1);
            pair<LL, ID> longest = { 0, 0 };
            auto fv = [&](ID s){
                if (covered[s] == 0)
                    dis[s] += VL[s].length();
                longest = max(longest, { dis[s], s });
            };
            auto fe = [&](ID s, ID t){
                if (dis[t] <= dis[s]) {
                    dis[t] = dis[s];
                    pre[t] = s;
                }
            };
            topo_sort_fn(fv, fe);
            vector<ID> path;
            ID s = longest.second;
            for (; pre[s] != -1; s = pre[s])
                path.push_back(s);
            path.push_back(s);
            reverse(path.begin(), path.end());
            paths.push_back(path);
            LL new_cover_cnt = 0;
            for (ID s : path) {
                if (covered[s] == 0)
                    new_cover_cnt++;
                covered[s]++;
            }
            cover_cnt += new_cover_cnt;
            cout << longest.first << " bps, " << new_cover_cnt << " in " << path.size() << " nodes path ending at node " << longest.second
                 << ", now " << cover_cnt << " / " << VL.size() << " nodes covered" << endl;
            if (limit != -1 && paths.size() >= limit)
            	break;
        }
        cout << paths.size() << " paths to cover" << (cover_cnt < VL.size() ? " (incompletely)" : "") << endl;
        return paths;
    }

    vector<ID> topo_sort_order() {
    	vector<ID> ids;
    	ids.reserve(VL.size());
        topo_sort_fn([&](const ID &s){ ids.push_back(s); }, [&](const ID &s, const ID &t){ });
        return ids;
    }

    vector<ID> get_dag_cut_points() {
    	vector<ID> ret;
        reset_degree();
        queue<ID> Q;
        vector<LL> touch(VL.size(), 0);
        LL cut_edges = 0; // number of visited edges between poped nodes and other nodes
        for (ID i = 0; i < VL.size(); i++) 
            if (D[i] == 0)
                Q.push(i);
        while (!Q.empty()) {
            ID s = Q.front();
            Q.pop();
            cut_edges -= touch[s];
            if (!E[s].empty() && cut_edges == 0 && touch[s] != 0)
            	ret.push_back(s);
            for (ID t : E[s]) {
            	D[t]--; touch[t]++; cut_edges++;
                if (D[t] == 0)
                    Q.push(t);
            }
        }
        return ret;
    }

    vector<string> original_id;
    void reset_original_id() {
    	original_id.resize(VL.size());
    	for (const pair<string, ID> &p : id_map)
    		original_id[p.second] = p.first;
    }

    Graph subgraph(const vector<ID> &ids) {
    	reset_original_id();
    	unordered_set<ID> sub(ids.begin(), ids.end());
    	Graph G;
    	for (ID id : ids)
    		G.add_node(original_id[id], VL[id]);
    	for (ID s : ids)
    		for (ID t : E[s])
    			if (sub.count(t))
    				G.add_edge(original_id[s], original_id[t]);
    	return G;
    }

    // output in dot language
    void visualize(const std::string &filename) {
    	reset_original_id();
    	std::ofstream fout(filename);
    	fout << "digraph graphname {\nnode [shape=record];\nrankdir=LR;" << endl;
    	for (ID i = 0; i < VL.size(); i++)
    		fout << i << " [label=\"" << original_id[i] << "\\l" << VL[i] << "\\l\"]" << endl;
    	for (ID i = 0; i < VL.size(); i++) 
    		if (!E[i].empty()) {
    			fout << i << ":e -> {";
    			for (LL j = 0; j < E[i].size(); j++)
    				fout << E[i][j] << ":w" << (j == E[i].size() - 1 ? "}" : ", ");
    			fout << endl;
    		}
    	fout << "}\n";
    }
};

template<typename T>
struct Range {
	unordered_set<T> s;
	T min_x, max_x, sum_x;
	LL count;
	vector<T> data;
	Range() {}
	Range(const vector<T> &vec) {
		for (const T &x : vec) add(x);
	}
	void add(const T &x) {
		if (s.empty() || x < min_x)
			min_x = x;
		if (s.empty() || x > max_x)
			max_x = x;
		if (s.empty())
			count = 0;
		if (s.empty())
			sum_x = x;
		else
			sum_x += x;
		s.insert(x);
		count++;
		data.push_back(x);
	}
	string to_string() const {
		string ret = "";
		ret += "size="+std::to_string(count);
		ret += " unique="+std::to_string(s.size());
		T avg_x = sum_x / count;
		ret += " min="+std::to_string(min_x);
		ret += " max="+std::to_string(max_x);
		ret += " sum="+std::to_string(sum_x);
		ret += " avg="+std::to_string(avg_x);
		return ret;
	}
	static vector<pair<T, T>> compressed(vector<T> &data) {
		vector<pair<T, T>> ret;
		sort(data.begin(), data.end());
		for (size_t i = 0; i < data.size(); i++) {
			if (i + 2 < data.size() && data[i] + data[i + 2] == data[i + 1] + data[i + 1]) {
				size_t j = i + 2;
				while (j + 1 < data.size() && data[j - 1] + data[j + 1] == data[j] + data[j])
					j++;
				ret.push_back({ data[i], data[j] });
				i = j;
			}
			else
				ret.push_back({ data[i], data[i] });
		}
		return ret;
	}
	vector<pair<T, T>> compressed() {
		return Range<T>::compressed(data);
	}
	static string compressed_str(vector<T> &data) {
		sort(data.begin(), data.end());
		string ret = "";
		for (const pair<T, T> &p : compressed(data)) {
			if (ret.length() > 0)
				ret += ", ";
			if (p.first == p.second)
				ret += std::to_string(p.first);
			else
				ret += "[" + std::to_string(p.first) + ".." + std::to_string(p.second) + "]";
		}
		return "{ " + ret + " }";
	}
	string compressed_str() {
		return Range<T>::compressed_str(data);
	}

};

#endif