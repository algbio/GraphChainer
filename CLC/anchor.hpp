#include <vector>
#include <fstream>
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