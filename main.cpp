#include <iostream>

#include <mpc/mpc.h>

using namespace lemon;


int main() {
    ListDigraph g;

    // Vertices
    ListDigraph::Node x0 = g.addNode();
    ListDigraph::Node x1 = g.addNode();
    ListDigraph::Node x2 = g.addNode();
    ListDigraph::Node x3 = g.addNode();
    ListDigraph::Node x4 = g.addNode();
    ListDigraph::Node x5 = g.addNode();
    ListDigraph::Node x6 = g.addNode();
    ListDigraph::Node x7 = g.addNode();

    // Edges
    g.addArc(x0, x1);
    g.addArc(x0, x2);
    g.addArc(x1, x3);
    g.addArc(x2, x3);
    g.addArc(x3, x4);
    g.addArc(x4, x5);
    g.addArc(x4, x6);
    g.addArc(x5, x7);
    g.addArc(x6, x7);


    std::cout << "MPC Min-Flow<MinCostFlow> reduction" << std::endl;
    for (auto path: MPC(g)) {
        for (auto vertex: path) {
            std::cout << g.id(vertex) << " , ";
        }
        std::cout << std::endl;
    }
    
    return 0;
}
