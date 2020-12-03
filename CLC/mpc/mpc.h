#ifndef GRAPHALGCHAINER_MPC_H
#define GRAPHALGCHAINER_MPC_H


#include <lemon/list_graph.h>


/*
 * Computes a minimum path cover by
 * reducing the problem to Min-flow<MinCostFlow>
 *
 * It assumes g is a DAG
 */
std::vector<std::vector<lemon::ListDigraph::Node>> MPC(lemon::ListDigraph& g);


#endif //GRAPHALGCHAINER_MPC_H
