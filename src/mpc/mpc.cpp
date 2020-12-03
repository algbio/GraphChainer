#include <mpc/mpc.h>

#include <lemon/cost_scaling.h>
#include <lemon/dfs.h>

using namespace lemon;



std::vector<std::vector<ListDigraph::Node>> MPC(ListDigraph& g) {

    // Build the Min-Flow network reduction
    ListDigraph red;

    ListDigraph::NodeMap<ListDigraph::Node> v_in(g);
    ListDigraph::NodeMap<ListDigraph::Node> v_out(g);
    ListDigraph::NodeMap<ListDigraph::Node> original(red);
    ListDigraph::ArcMap<int> cost(red, 0);
    ListDigraph::ArcMap<int> demand(red, 0);
    ListDigraph::NodeMap<int> supply(red, 0);

    ListDigraph::Node s = red.addNode();
    ListDigraph::Node t = red.addNode();
    ListDigraph::Arc st = red.addArc(s, t);
    supply[s] = countNodes(g);
    supply[t] = -countNodes(g);


    // Set split vertices
    for (ListDigraph::NodeIt v(g); v != INVALID; ++v) {
        v_in[v] = red.addNode();
        v_out[v] = red.addNode();
        original[v_in[v]] = v;
        original[v_out[v]] = v;

        ListDigraph::Arc split = red.addArc(v_in[v], v_out[v]);
        demand[split] = 1;

        ListDigraph::Arc sv = red.addArc(s, v_in[v]);
        cost[sv] = 1;

        red.addArc(v_out[v], t);
    }

    // Set edges connecting split vertices
    for (ListDigraph::ArcIt e(g); e != INVALID; ++e) {
        ListDigraph::Node u = g.source(e);
        ListDigraph::Node v = g.target(e);

        red.addArc(v_out[u], v_in[v]);
    }



    // Use CostScaling for solving the min-flow
    CostScaling<ListDigraph> cs(red);
    cs.lowerMap(demand).costMap(cost).supplyMap(supply).run();


    // Obtain Flow solution
    ListDigraph::ArcMap<int64_t> flowMap(red);
    cs.flowMap(flowMap);

    // Extract the Minimum Path Cover solution from the flow
    std::vector<std::vector<ListDigraph::Node>> path_cover;

    // Remove 0 flow edges and st
    std::vector<ListDigraph::Arc> removal_list = {st};
    for (ListDigraph::ArcIt e(red); e != INVALID; ++e) {
        if (flowMap[e] == 0) {
            removal_list.push_back(e);
        }
    }
    for (ListDigraph::Arc e : removal_list) {
        red.erase(e);
    }


    Dfs<ListDigraph> dfs(red);
    bool reachable = dfs.run(s, t);

    while (reachable) {
        std::vector<ListDigraph::Node> path;
        removal_list.clear();

        ListDigraph::Node v = t;
        ListDigraph::Arc e;
        while ((e = dfs.predArc(v)) != INVALID) {
            if (v != t && (path.empty() || path.back() != original[v]))
                path.push_back(original[v]);
            flowMap[e]--;
            if (flowMap[e] == 0) {
                removal_list.push_back(e);
            }
            v = red.source(e);
        }
        std::reverse(path.begin(), path.end());

        path_cover.push_back(path);
        for (ListDigraph::Arc e : removal_list) {
            red.erase(e);
        }

        dfs = Dfs<ListDigraph>(red);
        reachable = dfs.run(s, t);
    }

    return path_cover;
}