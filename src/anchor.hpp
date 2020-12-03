#include <vector>
#include <fstream>
#include <sstream>
#include <string>

typedef long long LL;

struct Anchor {
	std::vector<LL> path; // path on the graph, only start and end is used here
	LL x, y; // interval on the sequence [x, y]
};

typedef std::vector<Anchor> Chain;

LL coverage(const Chain &c) {
	LL ret = 0, last = -1;
	for (const Anchor &A : c) {
		ret += std::max(0LL, A.y - std::max(last, A.x - 1));
		last = std::max(last, A.y);
	}
	return ret;
}

void saveAnchors(const std::string &filename, const std::vector<Anchor> &anchors) {
    std::ofstream fout(filename);
    for (const Anchor &A : anchors) {
        fout << A.path.size() << " " << A.x << " " << A.y << std::endl;
        for (LL i = 0; i < A.path.size(); i++) {
            if (i && i % 80 == 0) fout << std::endl;
            fout << A.path[i] << " ";
        }
        fout << std::endl;
    }
}

std::vector<Anchor> loadAnchors(const std::string &filename) {
    std::vector<Anchor> ret;
    std::ifstream fin(filename);
    LL n;
    while (fin >> n) {
        Anchor A;
        fin >> A.x >> A.y;
        while (n--) {
            LL x;
            fin >> x;
            A.path.push_back(x);
        }
        ret.push_back(A);
    }
    return ret;
}


struct AnchorSet {
    std::string name;
    std::vector<Anchor> anchors;
};
// chains.txt from getchains.py
// format : >name l r
//          i_0 i_1 ... (negative if reversed)
std::vector<AnchorSet> loadAnchorSets(const std::string &filename) {
    std::vector<AnchorSet> ret;
    std::unordered_map<std::string, LL> names_ids_map;
    std::ifstream fin(filename);
    std::string s;
    while (std::getline(fin, s)) {
        if (s[0] == '>') {
            std::istringstream ssin(s);
            std::string name;
            Anchor A;
            ssin >> name >> A.x >> A.y;
            name = name.substr(1);
            if (!names_ids_map.count(name)) {
                names_ids_map[name] = ret.size();
                ret.resize(ret.size() + 1);
                ret.back().name = name;
            }
            if (std::getline(fin, s)) {
                std::istringstream ssin(s);
                LL x;
                bool rev = false;
                while (ssin >> x) {
                    if (x < 0)
                        rev = true;
                        // x = -x;
                    A.path.push_back(x);
                }
                if (!rev)
                    ret[names_ids_map[name]].anchors.push_back(A);
            }
        }
    }
    return ret;
}