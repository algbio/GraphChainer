#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <random>
#include <queue>
#include <unordered_set>
#include "AlignmentGraph.h"
#include "CommonUtils.h"
#include "ThreadReadAssertion.h"

AlignmentGraph dummy;

AlignmentGraph AlignmentGraph::DummyGraph()
{
	return dummy;
}

AlignmentGraph::AlignmentGraph() :
nodeLength(),
nodeLookup(),
nodeIDs(),
inNeighbors(),
nodeSequences(),
bpSize(0),
ambiguousNodeSequences(),
firstAmbiguous(std::numeric_limits<size_t>::max()),
DBGoverlap(0),
finalized(false)
{
}

size_t AlignmentGraph::getDBGoverlap() const
{
	return DBGoverlap;
}

void AlignmentGraph::ReserveNodes(size_t numNodes, size_t numSplitNodes)
{
	nodeSequences.reserve(numSplitNodes);
	ambiguousNodeSequences.reserve(numSplitNodes);
	nodeLookup.reserve(numNodes);
	nodeIDs.reserve(numSplitNodes);
	nodeLength.reserve(numSplitNodes);
	inNeighbors.reserve(numSplitNodes);
	outNeighbors.reserve(numSplitNodes);
	reverse.reserve(numSplitNodes);
	nodeOffset.reserve(numSplitNodes);
}

void AlignmentGraph::AddNode(int nodeId, const std::string& sequence, const std::string& name, bool reverseNode, const std::vector<size_t>& breakpoints)
{
	assert(firstAmbiguous == std::numeric_limits<size_t>::max());
	assert(!finalized);
	//subgraph extraction might produce different subgraphs with common nodes
	//don't add duplicate nodes
	if (nodeLookup.count(nodeId) != 0) return;
	originalNodeSize[nodeId] = sequence.size();
	originalNodeName[nodeId] = name;
	assert(breakpoints.size() >= 2);
	assert(breakpoints[0] == 0);
	assert(breakpoints.back() == sequence.size());
	for (size_t breakpoint = 1; breakpoint < breakpoints.size(); breakpoint++)
	{
		if (breakpoints[breakpoint] == breakpoints[breakpoint-1]) continue;
		assert(breakpoints[breakpoint] > breakpoints[breakpoint-1]);
		for (size_t offset = breakpoints[breakpoint-1]; offset < breakpoints[breakpoint]; offset += SPLIT_NODE_SIZE)
		{
			size_t size = SPLIT_NODE_SIZE;
			if (breakpoints[breakpoint] - offset < size) size = breakpoints[breakpoint] - offset;
			assert(size > 0);
			AddNode(nodeId, offset, sequence.substr(offset, size), reverseNode);
			if (offset > 0)
			{
				assert(outNeighbors.size() >= 2);
				assert(outNeighbors.size() == inNeighbors.size());
				assert(nodeIDs.size() == outNeighbors.size());
				assert(nodeOffset.size() == outNeighbors.size());
				assert(nodeIDs[outNeighbors.size()-2] == nodeIDs[outNeighbors.size()-1]);
				assert(nodeOffset[outNeighbors.size()-2] + nodeLength[outNeighbors.size()-2] == nodeOffset[outNeighbors.size()-1]);
				outNeighbors[outNeighbors.size()-2].push_back(outNeighbors.size()-1);
				inNeighbors[inNeighbors.size()-1].push_back(inNeighbors.size()-2);
			}
		}
	}
}

void AlignmentGraph::AddNode(int nodeId, int offset, const std::string& sequence, bool reverseNode)
{
	assert(firstAmbiguous == std::numeric_limits<size_t>::max());
	assert(!finalized);
	assert(sequence.size() <= SPLIT_NODE_SIZE);

	bpSize += sequence.size();
	nodeLookup[nodeId].push_back(nodeLength.size());
	nodeLength.push_back(sequence.size());
	nodeIDs.push_back(nodeId);
	inNeighbors.emplace_back();
	outNeighbors.emplace_back();
	reverse.push_back(reverseNode);
	nodeOffset.push_back(offset);
	NodeChunkSequence normalSeq;
	for (size_t i = 0; i < CHUNKS_IN_NODE; i++)
	{
		normalSeq[i] = 0;
	}
	AmbiguousChunkSequence ambiguousSeq;
	ambiguousSeq.A = 0;
	ambiguousSeq.C = 0;
	ambiguousSeq.G = 0;
	ambiguousSeq.T = 0;
	bool ambiguous = false;
	assert(sequence.size() <= sizeof(size_t)*8);
	for (size_t i = 0; i < sequence.size(); i++)
	{
		size_t chunk = i / BP_IN_CHUNK;
		assert(chunk < CHUNKS_IN_NODE);
		size_t offset = (i % BP_IN_CHUNK) * 2;
		switch(sequence[i])
		{
			case 'a':
			case 'A':
				ambiguousSeq.A |= ((size_t)1) << (i);
				normalSeq[chunk] |= ((size_t)0) << offset;
				break;
			case 'c':
			case 'C':
				ambiguousSeq.C |= ((size_t)1) << (i);
				normalSeq[chunk] |= ((size_t)1) << offset;
				break;
			case 'g':
			case 'G':
				ambiguousSeq.G |= ((size_t)1) << (i);
				normalSeq[chunk] |= ((size_t)2) << offset;
				break;
			case 't':
			case 'T':
			case 'u':
			case 'U':
				ambiguousSeq.T |= ((size_t)1) << (i);
				normalSeq[chunk] |= ((size_t)3) << offset;
				break;
			case 'r':
			case 'R':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'y':
			case 'Y':
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 's':
			case 'S':
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'w':
			case 'W':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'k':
			case 'K':
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'm':
			case 'M':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'b':
			case 'B':
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'd':
			case 'D':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'h':
			case 'H':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'v':
			case 'V':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			case 'n':
			case 'N':
				ambiguousSeq.A |= ((size_t)1) << (i);
				ambiguousSeq.C |= ((size_t)1) << (i);
				ambiguousSeq.G |= ((size_t)1) << (i);
				ambiguousSeq.T |= ((size_t)1) << (i);
				ambiguous = true;
				break;
			default:
				assert(false);
		}
	}
	ambiguousNodes.push_back(ambiguous);
	if (ambiguous)
	{
		ambiguousNodeSequences.emplace_back(ambiguousSeq);
	}
	else
	{
		nodeSequences.emplace_back(normalSeq);
	}
	assert(nodeIDs.size() == nodeLength.size());
	assert(nodeLength.size() == inNeighbors.size());
	assert(inNeighbors.size() == outNeighbors.size());
}

void AlignmentGraph::AddEdgeNodeId(int node_id_from, int node_id_to, size_t startOffset)
{
	assert(firstAmbiguous == std::numeric_limits<size_t>::max());
	assert(!finalized);
	assert(nodeLookup.count(node_id_from) > 0);
	assert(nodeLookup.count(node_id_to) > 0);
	size_t from = nodeLookup.at(node_id_from).back();
	size_t to = std::numeric_limits<size_t>::max();
	assert(nodeOffset[from] + nodeLength[from] == originalNodeSize[node_id_from]);
	for (auto node : nodeLookup[node_id_to])
	{
		if (nodeOffset[node] == startOffset)
		{
			to = node;
		}
	}
	assert(to != std::numeric_limits<size_t>::max());
	//don't add double edges
	if (std::find(inNeighbors[to].begin(), inNeighbors[to].end(), from) == inNeighbors[to].end()) inNeighbors[to].push_back(from);
	if (std::find(outNeighbors[from].begin(), outNeighbors[from].end(), to) == outNeighbors[from].end()) outNeighbors[from].push_back(to);
}

void AlignmentGraph::Finalize(int wordSize)
{
	assert(nodeSequences.size() + ambiguousNodeSequences.size() == nodeLength.size());
	assert(inNeighbors.size() == nodeLength.size());
	assert(outNeighbors.size() == nodeLength.size());
	assert(reverse.size() == nodeLength.size());
	assert(nodeIDs.size() == nodeLength.size());
	RenumberAmbiguousToEnd();
	ambiguousNodes.clear();
	findLinearizable();
	doComponentOrder();
	findChains();
	std::cout << nodeLookup.size() << " original nodes, " << (nodeLookup.size()/2) << " in one strand" << std::endl;
	std::cout << nodeLength.size() << " split nodes, " << (nodeLength.size()/2) << " in one strand" << std::endl;
	std::cout << ambiguousNodeSequences.size() << " ambiguous split nodes" << std::endl;
	finalized = true;
	int specialNodes = 0;
	size_t edges = 0;
	size_t bps = 0;
	for (size_t i = 0; i < inNeighbors.size(); i++)
	{
		inNeighbors[i].shrink_to_fit();
		outNeighbors[i].shrink_to_fit();
		if (inNeighbors[i].size() >= 2) specialNodes++;
		edges += inNeighbors[i].size();
		bps += nodeLength[i];
	}
	std::cout << bps << " base pairs (total characters), " << (bps/2) << " in one strand" << std::endl;
	std::cout << edges << " edges, " << (edges/2) << " in one strand" << std::endl;
	std::cout << specialNodes << " nodes with in-degree >= 2, " << (specialNodes/2) << " in one strand" << std::endl;
	assert(nodeSequences.size() + ambiguousNodeSequences.size() == nodeLength.size());
	assert(inNeighbors.size() == nodeLength.size());
	assert(outNeighbors.size() == nodeLength.size());
	assert(reverse.size() == nodeLength.size());
	assert(nodeIDs.size() == nodeLength.size());
	assert(nodeOffset.size() == nodeLength.size());
	nodeLength.shrink_to_fit();
	nodeIDs.shrink_to_fit();
	inNeighbors.shrink_to_fit();
	outNeighbors.shrink_to_fit();
	reverse.shrink_to_fit();
	nodeSequences.shrink_to_fit();
	ambiguousNodeSequences.shrink_to_fit();
#ifndef NDEBUG
	for (auto pair : nodeLookup)
	{
		for (size_t i = 1; i < pair.second.size(); i++)
		{
			assert(nodeOffset[pair.second[i-1]] < nodeOffset[pair.second[i]]);
		}
	}
#endif
}

std::pair<bool, size_t> AlignmentGraph::findBubble(const size_t start, const std::vector<bool>& ignorableTip)
{
	std::vector<size_t> S;
	S.push_back(start);
	std::unordered_set<size_t> visited;
	std::unordered_set<size_t> seen;
	seen.insert(start);
	while (S.size() > 0)
	{
		const size_t v = S.back();
		S.pop_back();
		assert(seen.count(v) == 1);
		seen.erase(v);
		assert(visited.count(v) == 0);
		visited.insert(v);
		if (outNeighbors[v].size() == 0) return std::make_pair(false, 0);
		for (const size_t u : outNeighbors[v])
		{
			if (ignorableTip[u]) continue;
			if (u == v) continue;
			if (u == start) return std::make_pair(false, 0);
			assert(visited.count(u) == 0);
			seen.insert(u);
			bool hasNonvisitedParent = false;
			for (const size_t w : inNeighbors[u])
			{
				if (w == u) continue;
				if (!ignorableTip[w] && visited.count(w) == 0)
				{
					hasNonvisitedParent = true;
					break;
				}
			}
			if (!hasNonvisitedParent) S.push_back(u);
		}
		if (S.size() == 1 && seen.size() == 1 && seen.count(S[0]) == 1)
		{
			const size_t t = S.back();
			for (const size_t u : outNeighbors[t])
			{
				if (u == start) return std::make_pair(false, 0);
			}
			return std::make_pair(true, t);
		}
	}
	return std::make_pair(false, 0);
}

size_t find(std::vector<size_t>& parent, size_t item)
{
	if (parent[item] == item) return item;
	std::vector<size_t> stack;
	stack.push_back(item);
	while (parent[stack.back()] != stack.back()) stack.push_back(parent[stack.back()]);
	for (size_t i : stack) parent[i] = stack.back();
	return stack.back();
}

void merge(std::vector<size_t>& parent, std::vector<size_t>& rank, size_t left, size_t right)
{
	left = find(parent, left);
	right = find(parent, right);
	if (rank[left] < rank[right])
	{
		std::swap(left, right);
	}
	parent[right] = left;
	if (rank[left] == rank[right]) rank[left] += 1;
}

void AlignmentGraph::chainBubble(const size_t start, const std::vector<bool>& ignorableTip, std::vector<size_t>& rank)
{
	bool hasBubble;
	size_t bubbleEnd;
	std::tie(hasBubble, bubbleEnd) = findBubble(start, ignorableTip);
	if (!hasBubble) return;
	std::unordered_set<size_t> visited;
	std::vector<size_t> stack;
	stack.push_back(start);
	visited.insert(start);
	merge(chainNumber, rank, start, bubbleEnd);
	while (stack.size() > 0)
	{
		const size_t top = stack.back();
		stack.pop_back();
		if (visited.count(top) == 1) continue;
		if (ignorableTip[top]) continue;
		visited.insert(top);
		merge(chainNumber, rank, start, top);
		for (const auto neighbor : outNeighbors[top])
		{
			if (visited.count(neighbor) == 1) continue;
			if (neighbor == bubbleEnd) continue;
			stack.push_back(neighbor);
		}
	}
}

void AlignmentGraph::fixChainApproxPos(const size_t start)
{
	assert(chainApproxPos[start] == std::numeric_limits<size_t>::max());
	assert(std::numeric_limits<size_t>::max() / SPLIT_NODE_SIZE > nodeLength.size());
	assert(std::numeric_limits<size_t>::max() > nodeLength.size() * SPLIT_NODE_SIZE * 2);
	std::vector<std::pair<size_t, size_t>> stack;
	size_t chain = chainNumber[start];
	stack.emplace_back(start, (nodeLength.size() + 5) * SPLIT_NODE_SIZE);
	while (stack.size() > 0)
	{
		size_t v;
		size_t dist;
		std::tie(v, dist) = stack.back();
		stack.pop_back();
		if (chainApproxPos[v] != std::numeric_limits<size_t>::max()) continue;
		chainApproxPos[v] = dist;
		for (const size_t u : outNeighbors[v])
		{
			if (chainNumber[u] != chain) continue;
			if (chainApproxPos[u] != std::numeric_limits<size_t>::max()) continue;
			assert(std::numeric_limits<size_t>::max() - nodeLength[u] > dist);
			stack.emplace_back(u, dist + nodeLength[u]);
		}
		for (const size_t u : inNeighbors[v])
		{
			if (chainNumber[u] != chain) continue;
			if (chainApproxPos[u] != std::numeric_limits<size_t>::max()) continue;
			assert(dist > nodeLength[v]);
			stack.emplace_back(u, dist - nodeLength[v]);
		}
	}
}

phmap::flat_hash_map<size_t, std::unordered_set<size_t>> AlignmentGraph::chainTips(std::vector<size_t>& rank, std::vector<bool>& ignorableTip)
{
	assert(componentNumber.size() == NodeSize());
	std::vector<size_t> order;
	order.reserve(NodeSize());
	for (size_t i = 0; i < NodeSize(); i++)
	{
		order.push_back(i);
	}
	std::sort(order.begin(), order.end(), [this](size_t left, size_t right) { return componentNumber[left] < componentNumber[right]; });
	std::vector<bool> fwTipComponent;
	fwTipComponent.resize(componentNumber[order.back()]+1, true);
	for (size_t ind = order.size()-1; ind < order.size(); ind--)
	{
		size_t i = order[ind];
		if (!fwTipComponent[componentNumber[i]]) continue;
		for (auto neighbor : outNeighbors[i])
		{
			assert(componentNumber[neighbor] >= componentNumber[i]);
			if (componentNumber[neighbor] == componentNumber[i])
			{
				fwTipComponent[componentNumber[i]] = false;
				break;
			}
			if (!fwTipComponent[componentNumber[neighbor]])
			{
				fwTipComponent[componentNumber[i]] = false;
				break;
			}
		}
	}
	for (size_t ind = order.size()-1; ind < order.size(); ind--)
	{
		size_t i = order[ind];
		if (!fwTipComponent[componentNumber[i]]) continue;
		for (auto neighbor : outNeighbors[i])
		{
			assert(fwTipComponent[componentNumber[neighbor]]);
			merge(chainNumber, rank, i, neighbor);
		}
	}
	std::vector<bool> bwTipComponent;
	bwTipComponent.resize(componentNumber[order.back()]+1, true);
	for (size_t ind = 0; ind < order.size(); ind++)
	{
		size_t i = order[ind];
		if (!bwTipComponent[componentNumber[i]]) continue;
		for (auto neighbor : inNeighbors[i])
		{
			assert(componentNumber[neighbor] <= componentNumber[i]);
			if (componentNumber[neighbor] == componentNumber[i])
			{
				bwTipComponent[componentNumber[i]] = false;
				break;
			}
			if (!bwTipComponent[componentNumber[neighbor]])
			{
				bwTipComponent[componentNumber[i]] = false;
				break;
			}
		}
	}
	for (size_t ind = 0; ind < order.size(); ind++)
	{
		size_t i = order[ind];
		if (!bwTipComponent[componentNumber[i]]) continue;
		for (auto neighbor : inNeighbors[i])
		{
			assert(bwTipComponent[componentNumber[neighbor]]);
			merge(chainNumber, rank, i, neighbor);
		}
	}
	phmap::flat_hash_map<size_t, std::unordered_set<size_t>> result;
	for (size_t i = 0; i < NodeSize(); i++)
	{
		if (bwTipComponent[componentNumber[i]] || fwTipComponent[componentNumber[i]])
		{
			ignorableTip[i] = true;
		}
		if (bwTipComponent[componentNumber[i]])
		{
			for (auto neighbor : outNeighbors[i])
			{
				if (chainNumber[neighbor] == chainNumber[i]) continue;
				result[chainNumber[i]].insert(neighbor);
			}
		}
		if (fwTipComponent[componentNumber[i]])
		{
			for (auto neighbor : inNeighbors[i])
			{
				if (chainNumber[neighbor] == chainNumber[i]) continue;
				result[chainNumber[i]].insert(neighbor);
			}
		}
	}
	return result;
}

void AlignmentGraph::chainCycles(std::vector<size_t>& rank, std::vector<bool>& ignorableTip)
{
	for (size_t i = 0; i < nodeLength.size(); i++)
	{
		size_t uniqueFwNeighbor = std::numeric_limits<size_t>::max();
		for (auto u : outNeighbors[i])
		{
			if (ignorableTip[u]) continue;
			if (u == i) continue;
			if (uniqueFwNeighbor == std::numeric_limits<size_t>::max())
			{
				uniqueFwNeighbor = u;
			}
			else
			{
				assert(u != uniqueFwNeighbor);
				uniqueFwNeighbor = std::numeric_limits<size_t>::max()-1;
			}
		}
		size_t uniqueBwNeighbor = std::numeric_limits<size_t>::max();
		for (auto u : inNeighbors[i])
		{
			if (ignorableTip[u]) continue;
			if (u == i) continue;
			if (uniqueBwNeighbor == std::numeric_limits<size_t>::max())
			{
				uniqueBwNeighbor = u;
			}
			else if (u != uniqueBwNeighbor)
			{
				uniqueBwNeighbor = std::numeric_limits<size_t>::max()-1;
			}
		}
		if (uniqueFwNeighbor != uniqueBwNeighbor) continue;
		if (uniqueFwNeighbor == std::numeric_limits<size_t>::max()) continue;
		if (uniqueFwNeighbor == std::numeric_limits<size_t>::max()-1) continue;
		if (uniqueBwNeighbor == std::numeric_limits<size_t>::max()) continue;
		if (uniqueBwNeighbor == std::numeric_limits<size_t>::max()-1) continue;
		ignorableTip[i] = true;
		assert(uniqueBwNeighbor == uniqueFwNeighbor);
		merge(chainNumber, rank, i, uniqueFwNeighbor);
	}
}

void AlignmentGraph::findChains()
{
	chainNumber.resize(nodeLength.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < chainNumber.size(); i++)
	{
		chainNumber[i] = i;
	}
	std::vector<bool> ignorableTip;
	ignorableTip.resize(nodeLength.size(), false);
	std::vector<size_t> rank;
	rank.resize(nodeLength.size(), 0);
	for (const auto& pair : nodeLookup)
	{
		assert(pair.second.size() > 0);
		for (size_t i = 1; i < pair.second.size(); i++)
		{
			merge(chainNumber, rank, pair.second[0], pair.second[i]);
		}
	}
	auto tipChainers = chainTips(rank, ignorableTip);
	chainCycles(rank, ignorableTip);
	for (const auto& pair : nodeLookup)
	{
		chainBubble(pair.second.back(), ignorableTip, rank);
	}
	for (auto& pair : tipChainers)
	{
		size_t uniqueNeighbor = std::numeric_limits<size_t>::max();
		for (auto n : pair.second)
		{
			if (uniqueNeighbor == std::numeric_limits<size_t>::max())
			{
				uniqueNeighbor = chainNumber[n];
			}
			if (uniqueNeighbor != chainNumber[n])
			{
				uniqueNeighbor = std::numeric_limits<size_t>::max()-1;
				break;
			}
		}
		assert(uniqueNeighbor != std::numeric_limits<size_t>::max());
		if (uniqueNeighbor == std::numeric_limits<size_t>::max()-1) continue;
		merge(chainNumber, rank, pair.first, *pair.second.begin());
	}
	{
		std::vector<size_t> tmp;
		std::vector<bool> tmp2;
		std::swap(rank, tmp);
		std::swap(ignorableTip, tmp2);
	}
	for (size_t i = 0; i < chainNumber.size(); i++)
	{
		find(chainNumber, i);
	}
	chainApproxPos.resize(nodeLength.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < chainNumber.size(); i++)
	{
		if (chainApproxPos[i] == std::numeric_limits<size_t>::max()) fixChainApproxPos(i);
	}
}

void AlignmentGraph::findLinearizable()
{
	linearizable.resize(nodeLength.size(), false);
	std::vector<bool> checked;
	checked.resize(nodeLength.size(), false);
	std::vector<size_t> stack;
	std::vector<bool> onStack;
	onStack.resize(nodeLength.size(), false);
	for (size_t node = 0; node < nodeLength.size(); node++)
	{
		if (checked[node]) continue;
		if (inNeighbors[node].size() != 1)
		{
			checked[node] = true;
			continue;
		}
		checked[node] = true;
		assert(inNeighbors[node].size() == 1);
		assert(stack.size() == 0);
		stack.push_back(node);
		onStack[node] = true;
		while (true)
		{
			assert(stack.size() <= nodeLength.size());
			if (inNeighbors[stack.back()].size() != 1)
			{
				for (size_t i = 0; i < stack.size()-1; i++)
				{
					assert(inNeighbors[stack[i]].size() == 1);
					checked[stack[i]] = true;
					linearizable[stack[i]] = true;
					onStack[stack[i]] = false;
				}
				linearizable[stack.back()] = false;
				checked[stack.back()] = true;
				onStack[stack.back()] = false;
				stack.clear();
				break;
			}
			assert(inNeighbors[stack.back()].size() == 1);
			if (checked[stack.back()])
			{
				for (size_t i = 0; i < stack.size()-1; i++)
				{
					assert(inNeighbors[stack[i]].size() == 1);
					checked[stack[i]] = true;
					linearizable[stack[i]] = true;
					onStack[stack[i]] = false;
				}
				linearizable[stack.back()] = false;
				checked[stack.back()] = true;
				onStack[stack.back()] = false;
				stack.clear();
				break;
			}
			assert(inNeighbors[stack.back()].size() == 1);
			auto neighbor = inNeighbors[stack.back()][0];
			if (neighbor == node)
			{
				for (size_t i = 0; i < stack.size(); i++)
				{
					checked[stack[i]] = true;
					linearizable[stack[i]] = false;
					onStack[stack[i]] = false;
				}
				stack.clear();
				break;
			}
			if (onStack[neighbor])
			{
				assert(neighbor != node);
				size_t i = stack.size();
				for (; i > 0; i--)
				{
					if (stack[i] == neighbor) break;
					checked[stack[i]] = true;
					linearizable[stack[i]] = false;
					onStack[stack[i]] = false;
				}
				for (size_t j = 0; j < i; j++)
				{
					checked[stack[j]] = true;
					linearizable[stack[j]] = true;
					onStack[stack[j]] = false;
				}
				stack.clear();
				break;
			}
			stack.push_back(inNeighbors[stack.back()][0]);
			onStack[stack.back()] = true;
		}
	}
}

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
size_t AlignmentGraph::NodeLength(size_t index) const
{
	return nodeLength[index];
}
size_t AlignmentGraph::NodeOffset(size_t index) const
{
	return nodeOffset[index];
}
size_t AlignmentGraph::NodeID(size_t index) const
{
	return nodeIDs[index];
}

char AlignmentGraph::NodeSequences(size_t node, size_t pos) const
{
	assert(pos < nodeLength[node]);
	if (node < firstAmbiguous)
	{
		assert(node < nodeSequences.size());
		size_t chunk = pos / BP_IN_CHUNK;
		size_t offset = (pos % BP_IN_CHUNK) * 2;
		return "ACGT"[(nodeSequences[node][chunk] >> offset) & 3];
	}
	else
	{
		assert(node >= firstAmbiguous);
		assert(node - firstAmbiguous < ambiguousNodeSequences.size());
		assert(pos < sizeof(size_t) * 8);
		bool A = (ambiguousNodeSequences[node - firstAmbiguous].A >> pos) & 1;
		bool C = (ambiguousNodeSequences[node - firstAmbiguous].C >> pos) & 1;
		bool G = (ambiguousNodeSequences[node - firstAmbiguous].G >> pos) & 1;
		bool T = (ambiguousNodeSequences[node - firstAmbiguous].T >> pos) & 1;
		assert(A + C + G + T >= 1);
		assert(A + C + G + T <= 4);
		if ( A && !C && !G && !T) return 'A';
		if (!A &&  C && !G && !T) return 'C';
		if (!A && !C &&  G && !T) return 'G';
		if (!A && !C && !G &&  T) return 'T';
		if ( A && !C &&  G && !T) return 'R';
		if (!A &&  C && !G &&  T) return 'Y';
		if (!A &&  C &&  G && !T) return 'S';
		if ( A && !C && !G &&  T) return 'W';
		if (!A && !C &&  G &&  T) return 'K';
		if ( A &&  C && !G && !T) return 'M';
		if (!A &&  C &&  G &&  T) return 'B';
		if ( A && !C &&  G &&  T) return 'D';
		if ( A &&  C && !G &&  T) return 'H';
		if ( A &&  C &&  G && !T) return 'V';
		if ( A &&  C &&  G &&  T) return 'N';
		assert(false);
		return 'N';
	}
}

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
AlignmentGraph::NodeChunkSequence AlignmentGraph::NodeChunks(size_t index) const
{
	assert(index < nodeSequences.size());
	return nodeSequences[index];
}

#ifdef NDEBUG
	__attribute__((always_inline))
#endif
AlignmentGraph::AmbiguousChunkSequence AlignmentGraph::AmbiguousNodeChunks(size_t index) const
{
	assert(index >= firstAmbiguous);
	assert(index - firstAmbiguous < ambiguousNodeSequences.size());
	return ambiguousNodeSequences[index - firstAmbiguous];
}

size_t AlignmentGraph::NodeSize() const
{
	return nodeLength.size();
}

class NodeWithDistance
{
public:
	NodeWithDistance(size_t node, bool start, size_t distance) : node(node), start(start), distance(distance) {};
	bool operator>(const NodeWithDistance& other) const
	{
		return distance > other.distance;
	}
	size_t node;
	bool start;
	size_t distance;
};

size_t AlignmentGraph::GetUnitigNode(int nodeId, size_t offset) const
{
	const auto& nodes = nodeLookup.at(nodeId);
	assert(nodes.size() > 0);
	//guess the index
	size_t index = nodes.size() * ((double)offset / (double)originalNodeSize.at(nodeId));
	if (index >= nodes.size()) index = nodes.size()-1;
	//go to the exact index
	while (index < nodes.size()-1 && (nodeOffset[nodes[index]] + NodeLength(nodes[index]) <= offset)) index++;
	while (index > 0 && (nodeOffset[nodes[index]] > offset)) index--;
	assert(index != nodes.size());
	size_t result = nodes[index];
	assert(nodeIDs[result] == nodeId);
	assert(nodeOffset[result] <= offset);
	assert(nodeOffset[result] + NodeLength(result) > offset);
	return result;
}

std::pair<int, size_t> AlignmentGraph::GetReversePosition(int nodeId, size_t offset) const
{
	assert(nodeLookup.count(nodeId) == 1);
	const auto& nodes = nodeLookup.at(nodeId);
	size_t originalSize = originalNodeSize.at(nodeId);
	assert(offset < originalSize);
	size_t newOffset = originalSize - offset - 1;
	assert(newOffset < originalSize);
	int reverseNodeId;
	if (nodeId % 2 == 0)
	{
		reverseNodeId = (nodeId / 2) * 2 + 1;
	}
	else
	{
		reverseNodeId = (nodeId / 2) * 2;
	}
	return std::make_pair(reverseNodeId, newOffset);
}

AlignmentGraph::MatrixPosition::MatrixPosition(size_t node, size_t nodeOffset, size_t seqPos) :
	node(node),
	nodeOffset(nodeOffset),
	seqPos(seqPos)
{
}

bool AlignmentGraph::MatrixPosition::operator==(const AlignmentGraph::MatrixPosition& other) const
{
	return node == other.node && nodeOffset == other.nodeOffset && seqPos == other.seqPos;
}

bool AlignmentGraph::MatrixPosition::operator!=(const AlignmentGraph::MatrixPosition& other) const
{
	return !(*this == other);
}

std::string AlignmentGraph::OriginalNodeName(int nodeId) const
{
	auto found = originalNodeName.find(nodeId);
	if (found == originalNodeName.end()) return "";
	return found->second;
}

std::vector<size_t> renumber(const std::vector<size_t>& vec, const std::vector<size_t>& renumbering)
{
	std::vector<size_t> result;
	result.reserve(vec.size());
	for (size_t i = 0; i < vec.size(); i++)
	{
		assert(vec[i] < renumbering.size());
		result.push_back(renumbering[vec[i]]);
	}
	return result;
}

template <typename T>
std::vector<T> reorder(const std::vector<T>& vec, const std::vector<size_t>& renumbering)
{
	assert(vec.size() == renumbering.size());
	std::vector<T> result;
	result.resize(vec.size());
	for (size_t i = 0; i < vec.size(); i++)
	{
		result[renumbering[i]] = vec[i];
	}
	return result;
}

void AlignmentGraph::RenumberAmbiguousToEnd()
{
	assert(nodeSequences.size() + ambiguousNodeSequences.size() == nodeLength.size());
	assert(inNeighbors.size() == nodeLength.size());
	assert(outNeighbors.size() == nodeLength.size());
	assert(reverse.size() == nodeLength.size());
	assert(nodeIDs.size() == nodeLength.size());
	assert(ambiguousNodes.size() == nodeLength.size());
	assert(firstAmbiguous == std::numeric_limits<size_t>::max());
	assert(!finalized);
	std::vector<size_t> renumbering;
	renumbering.reserve(ambiguousNodes.size());
	size_t nonAmbiguousCount = 0;
	size_t ambiguousCount = 0;
	for (size_t i = 0; i < ambiguousNodes.size(); i++)
	{
		if (!ambiguousNodes[i])
		{
			renumbering.push_back(nonAmbiguousCount);
			nonAmbiguousCount++;
		}
		else
		{
			assert(ambiguousCount < ambiguousNodes.size());
			assert(ambiguousNodes.size()-1-ambiguousCount >= nonAmbiguousCount);
			renumbering.push_back(ambiguousNodes.size()-1-ambiguousCount);
			ambiguousCount++;
		}
	}
	assert(renumbering.size() == ambiguousNodes.size());
	assert(nonAmbiguousCount + ambiguousCount == ambiguousNodes.size());
	assert(ambiguousCount == ambiguousNodeSequences.size());
	assert(nonAmbiguousCount == nodeSequences.size());
	firstAmbiguous = nonAmbiguousCount;

	if (ambiguousCount == 0) return;

	//the ambiguous nodes were added in the reverse order, reverse the sequence containers too
	std::reverse(ambiguousNodeSequences.begin(), ambiguousNodeSequences.end());

	nodeLength = reorder(nodeLength, renumbering);
	nodeOffset = reorder(nodeOffset, renumbering);
	nodeIDs = reorder(nodeIDs, renumbering);
	inNeighbors = reorder(inNeighbors, renumbering);
	outNeighbors = reorder(outNeighbors, renumbering);
	reverse = reorder(reverse, renumbering);
	for (auto& pair : nodeLookup)
	{
		pair.second = renumber(pair.second, renumbering);
	}
	assert(inNeighbors.size() == outNeighbors.size());
	for (size_t i = 0; i < inNeighbors.size(); i++)
	{
		inNeighbors[i] = renumber(inNeighbors[i], renumbering);
		outNeighbors[i] = renumber(outNeighbors[i], renumbering);
	}

#ifndef NDEBUG
	assert(inNeighbors.size() == outNeighbors.size());
	for (size_t i = 0; i < inNeighbors.size(); i++)
	{
		for (auto neighbor : inNeighbors[i])
		{
			assert(std::find(outNeighbors[neighbor].begin(), outNeighbors[neighbor].end(), i) != outNeighbors[neighbor].end());
		}
		for (auto neighbor : outNeighbors[i])
		{
			assert(std::find(inNeighbors[neighbor].begin(), inNeighbors[neighbor].end(), i) != inNeighbors[neighbor].end());
		}
	}
	for (auto pair : nodeLookup)
	{
		size_t foundSize = 0;
		std::set<size_t> offsets;
		size_t lastOffset = 0;
		for (auto node : pair.second)
		{
			assert(offsets.count(nodeOffset[node]) == 0);
			assert(offsets.size() == 0 || nodeOffset[node] > lastOffset);
			lastOffset = nodeOffset[node];
			offsets.insert(nodeOffset[node]);
			assert(nodeIDs[node] == pair.first);
			foundSize += nodeLength[node];
		}
		assert(foundSize == originalNodeSize[pair.first]);
	}
#endif
}

void AlignmentGraph::doComponentOrder()
{
	std::vector<std::tuple<size_t, int, size_t>> callStack;
	size_t i = 0;
	std::vector<size_t> index;
	std::vector<size_t> lowlink;
	std::vector<bool> onStack;
	std::vector<size_t> stack;
	index.resize(nodeLength.size(), std::numeric_limits<size_t>::max());
	lowlink.resize(nodeLength.size(), std::numeric_limits<size_t>::max());
	onStack.resize(nodeLength.size(), false);
	size_t checknode = 0;
	size_t nextComponent = 0;
	componentNumber.resize(nodeLength.size(), std::numeric_limits<size_t>::max());
	while (true)
	{
		if (callStack.size() == 0)
		{
			while (checknode < nodeLength.size() && index[checknode] != std::numeric_limits<size_t>::max())
			{
				checknode++;
			}
			if (checknode == nodeLength.size()) break;
			callStack.emplace_back(checknode, 0, 0);
			checknode++;
		}
		auto top = callStack.back();
		const size_t v = std::get<0>(top);
		int state = std::get<1>(top);
		size_t w;
		size_t neighborI = std::get<2>(top);
		callStack.pop_back();
		switch(state)
		{
			case 0:
				assert(index[v] == std::numeric_limits<size_t>::max());
				assert(lowlink[v] == std::numeric_limits<size_t>::max());
				assert(!onStack[v]);
				index[v] = i;
				lowlink[v] = i;
				i += 1;
				stack.push_back(v);
				onStack[v] = true;
				[[fallthrough]];
			startloop:
			case 1:
				if (neighborI >= outNeighbors[v].size()) goto endloop;
				assert(neighborI < outNeighbors[v].size());
				w = outNeighbors[v][neighborI];
				if (index[w] == std::numeric_limits<size_t>::max())
				{
					assert(lowlink[w] == std::numeric_limits<size_t>::max());
					assert(!onStack[w]);
					callStack.emplace_back(v, 2, neighborI);
					callStack.emplace_back(w, 0, 0);
					continue;
				}
				else if (onStack[w])
				{
					lowlink[v] = std::min(lowlink[v], index[w]);
					neighborI += 1;
					goto startloop;
				}
				else
				{
					neighborI += 1;
					goto startloop;
				}
			case 2:
				assert(neighborI < outNeighbors[v].size());
				w = outNeighbors[v][neighborI];
				assert(index[w] != std::numeric_limits<size_t>::max());
				assert(lowlink[w] != std::numeric_limits<size_t>::max());
				lowlink[v] = std::min(lowlink[v], lowlink[w]);
				neighborI++;
				goto startloop;
			endloop:
			case 3:
				if (lowlink[v] == index[v])
				{
					do
					{
						w = stack.back();
						stack.pop_back();
						onStack[w] = false;
						componentNumber[w] = nextComponent;
					} while (w != v);
					nextComponent++;
				}
		}
	}
	assert(stack.size() == 0);
	for (size_t i = 0; i < componentNumber.size(); i++)
	{
		assert(componentNumber[i] != std::numeric_limits<size_t>::max());
		assert(componentNumber[i] <= nextComponent-1);
		componentNumber[i] = nextComponent-1-componentNumber[i];
	}
#ifdef EXTRACORRECTNESSASSERTIONS
	for (size_t i = 0; i < nodeLength.size(); i++)
	{
		for (auto neighbor : outNeighbors[i])
		{
			assert(componentNumber[neighbor] >= componentNumber[i]);
		}
	}
#endif
}

size_t AlignmentGraph::ComponentSize() const
{
	return componentNumber.size();
}

size_t AlignmentGraph::SizeInBP() const
{
	return bpSize;
}

typedef long long LL;

struct flowGraph {
	LL N, M, S, T;
	std::vector<LL> f, p, t, c;

	flowGraph(LL NN) : N(NN+2) {
		init(N);
		S = NN;
		T = NN + 1;
	}

	void init(LL N) {
		f.clear();
		f.resize(N, 0);
		t.clear();
		t.resize(2);
		p = t;
		c = t;
	}

	void add_edge(LL i, LL j, LL cap) {
		// std::cerr << "add edge " << i << " " << j << "  " << cap << std::endl;
		p.push_back(j);
		t.push_back(f[i]);
		c.push_back(cap);
		f[i] = t.size() - 1;
	}
};

std::vector<std::vector<size_t>> AlignmentGraph::shrink(size_t cid, const std::vector<std::vector<size_t>> &pc) {
	// graph should be DAG
	const std::vector<size_t> &cids = component_ids[cid];
	size_t N = cids.size();

	std::vector<std::vector<size_t>> ret;
	LL K = pc.size(), inf = pc.size();
	std::vector<LL> covered(N, 0), starts(N, 0), ends(N, 0);
	std::map<std::pair<LL, LL>, LL> edge_covered;
	for (auto path : pc) {
		for (LL i = 0; i < path.size(); i++) {
			covered[component_idx[path[i]]]++;
			if (i > 0)
				edge_covered[{ component_idx[path[i - 1]], component_idx[path[i]] }]++;
		}
		starts[component_idx[path[0]]]++;
		ends[component_idx[path.back()]]++;
	}
	flowGraph fg(N * 2);
	// i_in = i, i_out = i + N
	// add r(i, j) = c(j,i) + f(i,j) - l(i,j)
	auto add = [&](LL i, LL j, LL cap, LL l, LL ff) {
		// std::cerr << "add edge " << i << " " << j << "  " << cap << "  " << l << " " <<ff << std::endl;
		fg.add_edge(i, j, 0 + ff - l);
		fg.add_edge(j, i, cap - ff);
	};
	for (LL i = 0; i < N; i++)
		for (size_t jid : outNeighbors[cids[i]]) {
			size_t j = component_idx[jid];
			LL ff = edge_covered.count({i, j}) ? edge_covered[{i, j}] : 0;
			add(i + N, j, inf, 0, ff);
		}
	for (LL i = 0; i < N; i++) {
		add(i, i + N, inf, 1, covered[i]);
		add(fg.S, i, inf, 0, starts[i]);
		add(i + N, fg.T, inf, 0, ends[i]);
	}
	
	LL total = inf;
	std::cerr << "start shrinkink from width " << total << std::endl;
	std::vector<LL> Q(fg.N, 0), pre(fg.N, -1), d(fg.N, 0);
	while (1) {
		LL Qsize = 0;
		Q[Qsize++] = fg.S;
		for (LL i = 0; i < fg.N; i++) {
			pre[i] = -1;
			d[i] = 0;
		}
		d[fg.S] = 1;
		for (LL idx = 0; idx < Qsize && d[fg.T] == 0;) {
			LL i = Q[idx++];
			for (LL e = fg.f[i]; e; e = fg.t[e]) {
				LL j = fg.p[e]; 
				if (fg.c[e] > 0 && d[j] == 0) {
					d[j] = 1;
					pre[j] = e;
					Q[Qsize++] = j;
				}
			}
		}
		if (d[fg.T] == 0) break;
		std::vector<size_t> tmp;
		LL flow = fg.c[pre[fg.T]];
		for (LL i = fg.T; ;) {
			tmp.push_back(i);
			LL e = pre[i];
			if (e == -1) break;
			flow = std::min(flow, fg.c[e]);
			i = fg.p[e ^ 1];
		}
		for (LL i = fg.T; ;) {
			LL e = pre[i];
			if (e == -1) break;
			fg.c[e] -= flow;
			fg.c[e^1] += flow;
			i = fg.p[e ^ 1];
		}
		if (flow == 0) exit(1);
		total -= flow;
		// std::cerr << "Now shrink by " << flow << " to " << total << std::endl;
	}
	// convert flow back to path cover
	// ret = pc;
	// ret.resize(total);
	for (LL itr = 0; itr < total; itr++) {
		std::vector<size_t> tmp;
		for (LL i = fg.S; i != fg.T; ) {
			if (0 <= i && i < N)
				tmp.push_back(cids[i]);
			LL nxt = -1;
			for (LL e = fg.f[i]; e; e = fg.t[e]) {
				LL j = fg.p[e];
				LL ff = fg.c[e] + ((i < N && i + N == j) ? 1 : 0);
				if ((e & 1) == 0 && ff > 0) {
					nxt = j;
					fg.c[e]--;
					break;
				}
			}
			if (nxt == -1) {
				std::cerr << i << " not found nxt" << std::endl;
				return ret;
			}
			i = nxt;
		}
		ret.push_back(tmp);
	}
	return ret;
}

std::vector<std::vector<size_t>> AlignmentGraph::greedyCover(size_t cid) const {
	const std::vector<size_t> &cids = component_ids[cid];
	size_t N = cids.size();

	std::vector<std::vector<size_t>> ret;
	std::vector<size_t> covered(N, 0);
	size_t covered_cnt = 0;
	std::vector<std::pair<size_t, size_t>> d(N);
	std::vector<size_t> incd(N), Q(N);
	while (covered_cnt < covered.size()) {
		size_t Qsize = 0;
		for (size_t i = 0; i < N; i++) {
			d[i] = std::make_pair(0, i);
			incd[i] = inNeighbors[cids[i]].size();
			if (incd[i] == 0)
				Q[Qsize++] = i;
		}
		std::pair<size_t, size_t> best = {0, 0};
		for (size_t i = 0; i < Qsize; ) {
			size_t s = Q[i++];
			if (covered[s] == 0)
				d[s].first++;
			best = std::max(best, {d[s].first, s});
			for (size_t tid : outNeighbors[cids[s]]) {
				size_t t = component_idx[tid];
				incd[t]--;
				d[t] = std::max(d[t], {d[s].first, s});
				if (incd[t] == 0)
					Q[Qsize++] = t;
			}
		}
		std::vector<size_t> tmp, path;
		for (size_t i = best.second; d[i].second != i || i != tmp.back(); i = d[i].second)
			tmp.push_back(i);
		std::reverse(tmp.begin(), tmp.end());
		size_t l = 0, r = tmp.size() - 1;
		while (covered[tmp[l]]) l++;
		while (covered[tmp[r]]) r--;
		size_t new_covered = 0;
		for (size_t i = l; i <= r; i++) {
			path.push_back(cids[tmp[i]]);
			if (covered[tmp[i]] == 0)
				new_covered++;
			covered[tmp[i]]++;
		}
		covered_cnt += new_covered;
		std::cout << "cid = " << cid << " path #" << mpc[cid].size() << " : " << path.size() << " " << new_covered << " " << (N - covered_cnt) << std::endl;
		ret.push_back(path);
	}
	return ret;
}

void AlignmentGraph::computeMPCIndex(size_t cid, const std::vector<std::vector<size_t>> &pc) {
	const std::vector<size_t> &cids = component_ids[cid];
	size_t N = cids.size();

	std::vector<std::vector<LL>> last2reach;
	LL K = pc.size();
	// forwards[cid].resize(N);
	backwards[cid].resize(N);
	last2reach.resize(N, std::vector<LL>(K, -1));
	paths[cid].resize(N);
	for (LL i = 0; i < K; i++)
		for (LL j = 0; j < pc[i].size(); j++) {
			size_t x = component_idx[pc[i][j]];
			last2reach[x][i] = j;
			paths[cid][x].push_back(i);
		}
	
	std::vector<LL> incd(N, 0), Q;
	for (LL i = 0; i < N; i++) {
		incd[i] = inNeighbors[cids[i]].size();
		if (incd[i] == 0)
			Q.push_back(i);
	}
	topo_ids[cid].clear();
	topo_ids[cid].resize(N);
	topo[cid].clear();
	for (LL i = 0; i < Q.size(); ) {
		LL s = Q[i++];
		for (size_t tid : outNeighbors[cids[s]]) {
			size_t t = component_idx[tid];
			incd[t]--;
			if (incd[t] == 0)
				Q.push_back(t);
		}
		topo_ids[cid][s] = topo[cid].size();
		topo[cid].push_back(s);
	}
	
	for (LL i : Q) {
		for (size_t jid : outNeighbors[cids[i]]) {
			size_t j = component_idx[jid];
			for (LL k = 0; k < K; k++)
				last2reach[j][k] = std::max(last2reach[j][k], last2reach[i][k]);
		}
	}
	for (LL i = 0; i < N; i++)
		for (LL k = 0; k < K; k++) {
			LL &idx = last2reach[i][k];
			if (idx != -1 && component_idx[pc[k][idx]] == i)
				idx--;
			
			if (idx != -1) {
				size_t x = component_idx[pc[k][idx]];
				// forwards[cid][idx].push_back({ i, k });
				backwards[cid][i].push_back({x, k});
			}
		}
	// for (LL i = 0; i < N; i++) {
	// 	if (backwards[i].size() > 1) {
	// 		std::sort(backwards[i].begin(), backwards[i].end());
	// 		backwards[i].erase(std::unqiue(backwards[i].begin(), backwards[i].end()), backwards[i].end());
	// 	}
	// }
}

bool AlignmentGraph::checkMinPathCover(const std::vector<std::vector<size_t>> &pc) {
	std::vector<LL> ids;
	for (LL i = 0; i < pc.size(); i++)
		ids.push_back(0);
	auto reachable = [&](LL S, LL T) {
		std::vector<LL> Q = { S }, d(NodeSize(), -1);
		d[S] = 1;
		for (LL x = 0; x < Q.size() && d[T] == -1; x++) {
			LL u = Q[x];
			// std::cerr << "? " << S << " "  << T << " : " << x << " : " << u << std::endl;
			for (size_t v : outNeighbors[u])
				if (d[v] == -1) {
					Q.push_back(v);
					d[v] = 1;
				}
		}
		return d[T] != -1;
	};
	while (1) {
		bool pushed = false;
		for (LL i = 0; i < pc.size(); i++) {
			for (LL j = 0; j < pc.size(); j++)
				while (i != j && ids[i] < pc[i].size() && reachable(pc[i][ids[i]], pc[j][ids[j]]))
					ids[i]++, pushed = true;
			if (ids[i] >= pc[i].size())
				return false;
		}
		if (!pushed)
			break;
	}
	for (LL i = 0; i < pc.size(); i++)
		for (LL j = 0; j < pc.size(); j++)
			if (i != j && reachable(pc[i][ids[i]], pc[j][ids[j]]))
				return false;
	return true;
}

void AlignmentGraph::buildComponentsMap() {
	component_map.clear();
	component_idx.clear();
	component_ids.clear();
	LL N = NodeSize();
	component_map.resize(N, N+1);
	component_idx.resize(N, N+1);
	std::vector<size_t> Q;
	for (size_t S = 0; S < N; S++) {
		if (component_map[S] != N + 1)
			continue;
		Q.clear();
		Q.push_back(S);
		size_t c = component_ids.size();
		component_map[S] = c;
		component_idx[S] = 0;
		for (size_t i = 0; i < Q.size(); ) {
			size_t s = Q[i++];
			for (size_t t : outNeighbors[s])
				if (component_map[t] == N + 1) {
					component_map[t] = c;
					component_idx[t] = Q.size();
					Q.push_back(t);
				}
			for (size_t t : inNeighbors[s])
				if (component_map[t] == N + 1) {
					component_map[t] = c;
					component_idx[t] = Q.size();
					Q.push_back(t);
				}
		}
		component_ids.push_back(Q);
	}
}

void AlignmentGraph::buildMPC() {
	std::cout << "Build MPC Index" << std::endl;
	buildComponentsMap();
	std::cout << component_ids.size() << " connected components" << std::endl;
	mpc.resize(component_ids.size());
	topo.resize(mpc.size());
	topo_ids.resize(mpc.size());
	paths.resize(mpc.size());
	// forwards.resize(mpc.size());
	backwards.resize(mpc.size());
	size_t tw = 0, mw = 0;
	for (size_t cid = 0; cid < component_ids.size(); cid++) {
		mpc[cid] = greedyCover(cid);
		std::cout << "cid = " << cid << " greedy width " << mpc[cid].size() << std::endl;
		mpc[cid] = shrink(cid, mpc[cid]);
		std::cout << "cid = " << cid << " optimal width " << mpc[cid].size() << std::endl;
		computeMPCIndex(cid, mpc[cid]);
		std::cout << "cid = " << cid << " MPC index done" << std::endl;
		tw += mpc[cid].size();
		mw = std::max(mw, mpc[cid].size());
	}
	std::cout << "MPC building done" << std::endl;
	std::cout << "total width " << tw << " and max component width " << mw << std::endl;
	// std::cout << checkMinPathCover(mpc) << std::endl;
}
void AlignmentGraph::loadMPC(const std::string &filename) {
	
}
void AlignmentGraph::saveMPC(const std::string &filename) {
	
}


std::vector<size_t> AlignmentGraph::generatePath(const std::string &seq_out, const std::string &path_out, const size_t seed) {
	// std::random_device rd;
	std::mt19937 gen(seed);
	gen.discard(700000); // discard first a few to properly initialize the random engine
	auto choose = [&](const std::vector<size_t> &vec) {
		return vec[gen() % vec.size()];
	};

	std::vector<size_t> path, Q;
	for (size_t i = 0; i < NodeSize(); i++)
		if (inNeighbors[i].size() == 0 && !reverse[i])
			Q.push_back(i);
	size_t s = choose(Q);
	while (1) {
		path.push_back(s);
		if (outNeighbors[s].empty())
			break;
		s = choose(outNeighbors[s]);
	}
	std::ofstream fout(seq_out), pout(path_out);
	fout << ">path_" << path[0] << "_" << path.back() << std::endl;
	int last_nodeID = -1; // Assumed not to be a valid ID
        for (size_t i : path) {
		if (nodeIDs[i] != last_nodeID) {
                	pout << i << ' ' << nodeIDs[i] << std::endl;
                       	last_nodeID = nodeIDs[i];
		}
		// std::cout << i << " " << nodeIDs[i] << " " << reverse[i] << " : ";
		for (size_t j = 0; j < NodeLength(i); j++) {
			fout << NodeSequences(i, j);
			// std::cout << NodeSequences(i, j);
		}
		// std::cout << std::endl;
	}
	return path;
}


template<typename T>
struct SegmentTree {
	typedef long long LL;
	LL N;
	std::vector<T> t;
	LL ql, qr;
	T qx;

	LL get_upper_bit(LL N) {
		LL x = N & -N;
		while (x != N) {
			N ^= x;
			x = N & -N;
		}
		return x << 1;
	}
	static LL ID(LL l, LL r) {
		return l + r - 1;
	}
	SegmentTree(LL N_, T value) : N(get_upper_bit(N_)) {
		reset(value);
	}
	void reset(T value) {
		t.clear();
		t.resize(N*2, value);
	}
	void add(LL x, T value) {
		ql = x; qx = value;
		add_(0, N);
	}
	void add_(LL l, LL r) {
		t[ID(l, r)] = std::max(t[ID(l, r)], qx);
		if (l + 1 == r)
			return;
		LL m = (l + r) >> 1;
		if (ql < m)
			add_(l, m);
		else
			add_(m, r);
	}
	T RMQ(LL l, LL r) {
		if (l > r)
			return t[0];
		ql = l; qr = r + 1;
		return RMQ_(0, N);
	}
	T RMQ_(LL l, LL r) {
		// if (ql==27579&&qr==27605+1) cout<<"now "<<ID(l,r)<<"  "<<l<<" "<<r<<" : "<<t[ID(l,r)].first<<" "<<t[ID(l,r)].second<<endl;
		if (ql <= l && r <= qr)
			return t[ID(l, r)];
		LL m = (l + r) >> 1;
		if (qr <= m)
			return RMQ_(l, m);
		if (m <= ql)
			return RMQ_(m, r);
		return std::max(RMQ_(l, m), RMQ_(m, r));
	}
};


template<typename T, typename V>
struct Treap {
	struct Node {
		int ls, rs, size, pri;
		T key;
		V value, max;
	};
	std::vector<Node> t;
	int root;
	V default_value;
	Treap(const V &default_value = V()) : default_value(default_value) {
		root = 0;
		t.resize(1);
	}
	inline int randomm() {
		static int seed = 703; 
		return seed = int(seed * 48271LL % 2147483647);
	}
	inline int update(int now) {
		t[now].size = 1;
		t[now].max = t[now].value;
		if (t[now].ls) {
			t[now].size += t[t[now].ls].size;
			t[now].max = max(t[now].max, t[t[now].ls].max);
		}
		if (t[now].rs) {
			t[now].size += t[t[now].rs].size;
			t[now].max = max(t[now].max, t[t[now].rs].max);
		}
		return now;
	}
	inline int new_node (T key, V value) {
		t.push_back(Node({ 0, 0, 1, randomm(), key, value, value }));
		return t.size() - 1;
	}
	int merge(int x, int y) {
		if (!x || !y) return x + y;
		if (t[x].pri > t[y].pri) {
			t[x].rs = merge(t[x].rs, y);
			return update(x);
		}
		else {
			t[y].ls = merge(x, t[y].ls);
			return update(y);
		}
	}
	void split(int now, T key, int &x, int &y) {
		if (!now) {
			x = y = 0; 
			return;
		}
		if (t[now].key <= key) {
			x = now;
			split(t[now].rs, key, t[now].rs, y);
			update(x);
		}
		else {
			y = now;
			split(t[now].ls, key, x, t[now].ls);
			update(y);
		}
	}
	// void Del(int &root, int key) {
	// 	int x = 0, y = 0, z = 0;
	// 	split(root, key, x, z);
	// 	split(x, key - 1, x, y);
	// 	y = merge(t[y].ls, t[y].rs);
	// 	root = merge(merge(x, y), z);
	// }
	void add(T key, V value) {
		int x = 0, y = 0, z = 0;
		split(root, key, x, y);
		root = merge(merge(x, new_node(key, value)), y);
	}	
	V RMQ(T l, T r) {
		int now = root;
		while (now != 0 && (t[now].key < l || t[now].key > r)) {
			if (t[now].key < l)
				now = t[now].rs;
			else
				now = t[now].ls;
		}
		if (now == 0) {
			return default_value;
		}
		V ret = t[now].value;
		int x = t[now].ls;
		while (x != 0) {
			if (t[x].key >= l) {
				ret = max(ret, t[x].value);
				if (t[x].rs != 0)
					ret = max(ret, t[t[x].rs].max);
				x = t[x].ls;
			}
			else
				x = t[x].rs;
		}
		int y = t[now].rs;
		while (y != 0) {
			if (t[y].key <= r) {
				ret = max(ret, t[y].value);
				if (t[y].ls != 0)
					ret = max(ret, t[t[y].ls].max);
				y = t[y].rs;
			}
			else
				y = t[y].ls;
		}
		return ret;
	}
};

std::vector<size_t> AlignmentGraph::colinearChaining(const std::vector<Anchor> &A, long long sep_limit) const {
	std::vector<std::pair<size_t, size_t>> cs(A.size());
	for (size_t i = 0; i < A.size(); i++) {
		cs[i].first = component_map[A[i].path.back()];
		cs[i].second = i;
	}
	std::sort(cs.begin(), cs.end());
	std::vector<size_t> aids;
	std::pair<std::vector<size_t>, size_t> best, tmp;
	bool first = true;
	for (size_t i = 0, j; i < cs.size(); ) {
		aids.clear();
		for (j = i; j < cs.size() && cs[j].first == cs[i].first; j++)
			aids.push_back(cs[j].second);
		tmp = colinearChainingByComponent(cs[i].first, A, aids, sep_limit);
		// std::cerr << "cid " << cs[i].first << " " << aids.size() << " / " << A.size() << " : " << tmp.second << std::endl;
		if (first || tmp.second > best.second) {
			first = false;
			best = tmp;
		}
		i = j;
	}
	return best.first;
}

std::pair<std::vector<size_t>, size_t> AlignmentGraph::colinearChainingByComponent(size_t cid, const std::vector<Anchor> &A, const std::vector<size_t> &aids, long long sep_limit) const {
	typedef long long LL;
	auto getSortedMap = [&](std::vector<LL> a) {
		std::sort(a.begin(), a.end());
		a.erase(std::unique(a.begin(), a.end()), a.end());
		std::unordered_map<LL, LL> ret;
		for (LL i = 0; i < a.size(); i++)
			ret[a[i]] = i;
		return ret;
	};
	const std::vector<size_t> &cids = component_ids[cid];
	size_t N = cids.size();

	LL K = mpc[cid].size();
	std::pair<LL, LL> defaul_value = { -N*2, -1 };
	for (size_t j : aids) {
		defaul_value.first -= (A[j].y + 1 - A[j].x) * 2;
	}
	// std::cerr <<"defaul_value "<<defaul_value.first<<std::endl;
	typedef Treap<LL, std::pair<LL, LL>> IndexT;
	// std::vector<SegmentTree<std::pair<LL, LL>>> T(K, SegmentTree(N, defaul_value)), I(K, SegmentTree(N, defaul_value));
	std::vector<IndexT> T(K, IndexT(defaul_value)), I(K, IndexT(defaul_value));
	// std::vector<std::vector<LL>> starts(N), ends(N);
	std::vector<std::pair<LL, std::pair<LL, LL>>> endpoints;
	std::vector<std::pair<LL, LL>> C(A.size());
	for (size_t j : aids) {
		endpoints.push_back({ component_idx[A[j].path[0]], {j, -1} });
		endpoints.push_back({ component_idx[A[j].path.back()], {j, -2} });
		for (std::pair<LL, LL> b : backwards[cid][component_idx[A[j].path[0]]]) {
			endpoints.push_back({b.first, {j, b.second}});
			// std::cerr << "  A " << j << " back " << b.first << " " << b.second << std::endl;
		}
		C[j] = { A[j].y - A[j].x + 1, -1 };
		// std::cerr << "  A " << j << " " << A[j].path[0] << " " << A[j].path[1] << " -> " << component_idx[A[j].path[0]] << " " << component_idx[A[j].path.back()] << std::endl;
	}
	std::sort(endpoints.begin(), endpoints.end(), [&](const std::pair<LL, std::pair<LL, LL>> &p1, const std::pair<LL, std::pair<LL, LL>> &p2){
		return topo_ids[cid][p1.first] < topo_ids[cid][p2.first];
	});
	// for (LL vidx = 0; vidx < N; vidx++) {
	// 	LL v = topo[vidx];
	for (LL vidx = 0, ridx = 0; vidx < endpoints.size(); vidx = ridx) {
		// if (vidx > 0 && endpoints[vidx - 1] == endpoints[vidx])
		// 	continue;
		LL v = endpoints[vidx].first;
		ridx = vidx + 1;
		while (ridx < endpoints.size() && endpoints[ridx].first == v)
			ridx++;
		// std::cerr <<"now v="<<v<< " topoidx "<< topo_ids[cid][v] << " " << vidx << " " << ridx << "  " << std::endl;
		std::vector<LL> ids;
		for (size_t j = vidx; j < ridx; j++)
			if (endpoints[j].second.second < 0)
				ids.push_back(endpoints[j].second.first);
		if (ids.size() > 0) {
			std::sort(ids.begin(), ids.end(), [&](LL i, LL j) { 
				if (A[i].y != A[j].y)
					return A[i].y < A[j].y;
				return A[i].x < A[j].x;
			});
			ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
			// get local id count
			std::vector<LL> pos = { 0 };
			for (LL j : ids) {
				pos.push_back(A[j].x - 1);
				pos.push_back(A[j].x);
				pos.push_back(A[j].y - 1);
				pos.push_back(A[j].y);
			}
			auto id_map = getSortedMap(pos);
			LL Size = (LL)id_map.size();
			// IndexT tmpT(Size, defaul_value), tmpI(Size, defaul_value);
			IndexT tmpT(defaul_value), tmpI(defaul_value);
			for (LL j : ids) {
				if (component_idx[A[j].path[0]] == v) {
					std::pair<LL, LL> q = tmpT.RMQ(id_map[0], id_map[A[j].x - 1]);
					// if (q.second!=-1)std::cerr << "C " << j << "updates A " << C[j].first << " " << C[j].second << " <- " << A[j].y - A[j].x + 1 + q.first << " " << q.second << std::endl;
					C[j] = std::max(C[j], {A[j].y - A[j].x + 1 + q.first, q.second});
					q = tmpI.RMQ(id_map[A[j].x], id_map[A[j].y - 1]);
					// if (q.second!=-1)std::cerr << "C " << j << "updates B " << C[j].first << " " << C[j].second << " <- " << A[j].y + q.first << " " << q.second << std::endl;
					C[j] = std::max(C[j], {A[j].y + q.first, q.second});
				}
				if (component_idx[A[j].path.back()] == v) {
					tmpT.add(id_map[A[j].y], {C[j].first, j});
					tmpI.add(id_map[A[j].y], {C[j].first - A[j].y, j});
				}
			}
		}
		for (size_t vi = vidx; vi < ridx; vi++) {
			if (endpoints[vi].second.second != -2)
				continue;
			size_t j = endpoints[vi].second.first;
			// if (v == component_idx[A[j].path.back()])std::cerr << "add to " << j << std::endl;
			if (v == component_idx[A[j].path.back()])
				for (LL k : paths[cid][v]) {
					T[k].add(A[j].y, {C[j].first, j});
					I[k].add(A[j].y, {C[j].first - A[j].y, j});
				}
		}
		for (size_t vi = vidx; vi < ridx; vi++) {
			if (endpoints[vi].second.second < 0)
				continue;
			size_t j = endpoints[vi].second.first;
			LL u = component_idx[A[j].path[0]], k = endpoints[vi].second.second;
			std::pair<LL, LL> q = T[k].RMQ(0, A[j].x - 1);
			// if (q.second!=-1)std::cerr << "C " << j << "updates C " << C[j].first << " " << C[j].second << " <- " << A[j].y - A[j].x + 1 + q.first << " " << q.second << std::endl;
			C[j] = std::max(C[j], {A[j].y - A[j].x + 1 + q.first, q.second}); 
			q = I[k].RMQ(A[j].x, A[j].y - 1);
			// if (q.second!=-1)std::cerr << "C " << j << "updates D " << C[j].first << " " << C[j].second << " <- " << A[j].y + q.first << " " << q.second << std::endl;
			C[j] = std::max(C[j], {A[j].y + q.first, q.second});
		}
	}
	std::pair<LL, LL> best = {0, -1};
	for (size_t j : aids) 
		best = std::max(best, { C[j].first, j });
	// std::cerr << "optimal coverage : " << best.first << std::endl;
	// std::cerr << "optimal ends : " << best.second << std::endl;
	std::vector<size_t> ret;
	for (LL i = best.second; i != -1; i = C[i].second) {
		ret.push_back(i);
		// std::cerr << "now " << i << "  " << C[i].first << " "  << C[i].second << endl;
		if (i == C[i].second) {
			std::cerr << "error, loops in C[j] : " << i << "  " << C[i].first << " "  << C[i].second << std::endl;
			break;
		}
	}
	std::reverse(ret.begin(), ret.end());
	return { ret, best.first };
}


std::vector<size_t> AlignmentGraph::getChainPath(size_t S, size_t T, long long sep_limit) const {
// std::vector<std::pair<size_t, std::pair<size_t, size_t>>> AlignmentGraph::getChainPath(const std::vector<Anchor> &anchors, const std::vector<size_t> &ids, long long sep_limit) const {
	// std::unordered_set<LL> nodes;
	// std::vector<std::pair<LL, std::pair<LL, std::pair<LL, LL>>>> ids;
	// for (size_t i : ids) {
	// 	for (size_t j = 0; j < anchors[i].path.size(); j++) {
	// 		size_t x = anchors[i].path[j];
	// 		if (!nodes.count(x)) {
	// 			ids.push_back({ topo_ids[x], {x, {i, j}} });
	// 			nodes.insert(x);
	// 		}
	// }
	// std::sort(ids.begin(), ids.end());
	// ids.erase( std::unique( ids.begin(), ids.end() ), ids.end() );
	// LL last = -1;
	LL N = NodeSize();
	// std::vector<std::pair<size_t, std::pair<size_t, size_t>>> ret, longest, tmp;
	thread_local std::vector<size_t> vis, dis, Q, pre, tmp;
	thread_local size_t flag = 1;
	if (vis.size() < N) {
		vis.resize(N, 0);
		pre.resize(N);
		dis.resize(N);
		Q.reserve(N);
	}
	Q.clear();
	Q.push_back(S);
	vis[S] = ++flag;
	dis[S] = 0;
	for (size_t i = 0; vis[T] != flag && i < Q.size(); ) {
		size_t s = Q[i++];
		if (dis[s] > sep_limit) {
			continue;
		}
		for (size_t t : outNeighbors[s])
			if (vis[t] != flag) {
				Q.push_back(t);
				vis[t] = flag;
				dis[t] = dis[s] + NodeLength(t);
				pre[t] = s;
			}
	}
	tmp.clear();
	if (vis[T] != flag)
		return tmp;
	for (size_t i = T; i != S; i = pre[i])
		tmp.push_back(i);
	tmp.push_back(S);
	std::reverse(tmp.begin(), tmp.end());
	return tmp;
}
