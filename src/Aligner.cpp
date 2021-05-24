#include <sstream>
#include <iostream>
#include <fstream>
#include <functional>
#include <algorithm>
#include <thread>
#include <concurrentqueue.h> //https://github.com/cameron314/concurrentqueue
#include <google/protobuf/util/json_util.h>
#include "Aligner.h"
#include "CommonUtils.h"
#include "vg.pb.h"
#include "stream.hpp"
#include "fastqloader.h"
#include "BigraphToDigraph.h"
#include "ThreadReadAssertion.h"
#include "GraphAlignerWrapper.h"
#include "MummerSeeder.h"
#include "ReadCorrection.h"
#include "MinimizerSeeder.h"
#include "AlignmentSelection.h"

#include <edlib.h>

struct Seeder
{
	enum Mode
	{
		None, File, Mum, Mem, Minimizer
	};
	Mode mode;
	size_t mumCount;
	size_t memCount;
	size_t mxmLength;
	size_t minimizerLength;
	size_t minimizerWindowSize;
	double minimizerSeedDensity;
	const MummerSeeder* mummerSeeder;
	const MinimizerSeeder* minimizerSeeder;
	const std::unordered_map<std::string, std::vector<SeedHit>>* fileSeeds;
	Seeder(const AlignerParams& params, const std::unordered_map<std::string, std::vector<SeedHit>>* fileSeeds, const MummerSeeder* mummerSeeder, const MinimizerSeeder* minimizerSeeder) :
		mumCount(params.mumCount),
		memCount(params.memCount),
		mxmLength(params.mxmLength),
		minimizerLength(params.minimizerLength),
		minimizerWindowSize(params.minimizerWindowSize),
		minimizerSeedDensity(params.minimizerSeedDensity),
		mummerSeeder(mummerSeeder),
		minimizerSeeder(minimizerSeeder),
		fileSeeds(fileSeeds)
	{
		mode = Mode::None;
		if (fileSeeds != nullptr)
		{
			assert(minimizerSeeder == nullptr);
			assert(mummerSeeder == nullptr);
			assert(mumCount == 0);
			assert(memCount == 0);
			assert(minimizerSeedDensity == 0);
			mode = Mode::File;
		}
		if (minimizerSeeder != nullptr)
		{
			assert(mummerSeeder == nullptr);
			assert(mumCount == 0);
			assert(memCount == 0);
			assert(minimizerSeedDensity != 0);
			mode = Mode::Minimizer;
		}
		if (mummerSeeder != nullptr)
		{
			assert(minimizerSeeder == nullptr);
			assert(fileSeeds == nullptr);
			assert(mumCount != 0 || memCount != 0);
			assert(minimizerSeedDensity == 0);
			if (mumCount != 0)
			{
				mode = Mode::Mum;
				assert(memCount == 0);
			}
			if (memCount != 0)
			{
				mode = Mode::Mem;
				assert(mumCount == 0);
			}
		}
	}
	std::vector<SeedHit> getSeeds(const std::string& seqName, const std::string& seq) const
	{
		switch(mode)
		{
			case Mode::File:
				assert(fileSeeds != nullptr);
				if (fileSeeds->count(seqName) == 0) return std::vector<SeedHit>{};
				return fileSeeds->at(seqName);
			case Mode::Mum:
				assert(mummerSeeder != nullptr);
				return mummerSeeder->getMumSeeds(seq, mumCount, mxmLength);
			case Mode::Mem:
				assert(mummerSeeder != nullptr);
				return mummerSeeder->getMemSeeds(seq, memCount, mxmLength);
			case Mode::Minimizer:
				assert(minimizerSeeder != nullptr);
				return minimizerSeeder->getSeeds(seq, minimizerSeedDensity);
			case Mode::None:
				assert(false);
		}
		return std::vector<SeedHit>{};
	}
};

struct AlignmentStats
{
	AlignmentStats() :
	reads(0),
	seeds(0),
	seedsFound(0),
	seedsExtended(0),
	readsWithASeed(0),
	alignments(0),
	fullLengthAlignments(0),
	readsWithAnAlignment(0),
	bpInReads(0),
	bpInReadsWithASeed(0),
	bpInAlignments(0),
	bpInFullAlignments(0),
	allAlignmentsCount(0),
	assertionBroke(false)
	{
	}
	std::atomic<size_t> reads;
	std::atomic<size_t> seeds;
	std::atomic<size_t> seedsFound;
	std::atomic<size_t> seedsExtended;
	std::atomic<size_t> readsWithASeed;
	std::atomic<size_t> alignments;
	std::atomic<size_t> fullLengthAlignments;
	std::atomic<size_t> readsWithAnAlignment;
	std::atomic<size_t> bpInReads;
	std::atomic<size_t> bpInReadsWithASeed;
	std::atomic<size_t> bpInAlignments;
	std::atomic<size_t> bpInFullAlignments;
	std::atomic<size_t> allAlignmentsCount;
	std::atomic<bool> assertionBroke;
};

bool is_file_exist(std::string fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}

void replaceDigraphNodeIdsWithOriginalNodeIds(vg::Alignment& alignment, const AlignmentGraph& graph)
{
	for (int i = 0; i < alignment.path().mapping_size(); i++)
	{
		int digraphNodeId = alignment.path().mapping(i).position().node_id();
		int originalNodeId = digraphNodeId / 2;
		alignment.mutable_path()->mutable_mapping(i)->mutable_position()->set_node_id(originalNodeId);
		std::string name = graph.OriginalNodeName(digraphNodeId);
		if (name.size() > 0)
		{
			alignment.mutable_path()->mutable_mapping(i)->mutable_position()->set_name(name);
		}
	}
}

void readFastqs(const std::vector<std::string>& filenames, moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>>& writequeue, std::atomic<bool>& readStreamingFinished)
{
	assertSetNoRead("Read streamer");
	for (auto filename : filenames)
	{
		FastQ::streamFastqFromFile(filename, false, [&writequeue](FastQ& read)
		{
			std::shared_ptr<FastQ> ptr = std::make_shared<FastQ>();
			std::swap(*ptr, read);
			size_t slept = 0;
			while (writequeue.size_approx() > 200)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				slept++;
				if (slept > 100) break;
			}
			writequeue.enqueue(ptr);
		});
	}
	readStreamingFinished = true;
}

void consumeBytesAndWrite(const std::string& filename, moodycamel::ConcurrentQueue<std::string*>& writequeue, moodycamel::ConcurrentQueue<std::string*>& deallocqueue, std::atomic<bool>& allThreadsDone, std::atomic<bool>& allWriteDone, bool verboseMode, bool textMode)
{
	assertSetNoRead("Writer");
	auto openmode = std::ios::out;
	if (!textMode) openmode |= std::ios::binary;
	std::ofstream outfile { filename, openmode };

	bool wroteAny = false;

	std::string* alns[100] {};

	BufferedWriter coutoutput;
	if (verboseMode)
	{
		coutoutput = {std::cout};
	}

	while (true)
	{
		size_t gotAlns = writequeue.try_dequeue_bulk(alns, 100);
		if (gotAlns == 0)
		{
			if (!writequeue.try_dequeue(alns[0]))
			{
				if (allThreadsDone) break;
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				continue;
			}
			gotAlns = 1;
		}
		coutoutput << "write " << gotAlns << ", " << writequeue.size_approx() << " left" << BufferedWriter::Flush;
		for (size_t i = 0; i < gotAlns; i++)
		{
			outfile.write(alns[i]->data(), alns[i]->size());
		}
		deallocqueue.enqueue_bulk(alns, gotAlns);
		wroteAny = true;
	}

	if (!textMode && !wroteAny)
	{
		::google::protobuf::io::ZeroCopyOutputStream *raw_out =
		      new ::google::protobuf::io::OstreamOutputStream(&outfile);
		::google::protobuf::io::GzipOutputStream *gzip_out =
		      new ::google::protobuf::io::GzipOutputStream(raw_out);
		::google::protobuf::io::CodedOutputStream *coded_out =
		      new ::google::protobuf::io::CodedOutputStream(gzip_out);
		coded_out->WriteVarint64(0);
		delete coded_out;
		delete gzip_out;
		delete raw_out;
	}

	allWriteDone = true;
}

void QueueInsertSlowly(moodycamel::ProducerToken& token, moodycamel::ConcurrentQueue<std::string*>& queue, std::string&& str)
{
	std::string* write = new std::string { std::move(str) };
	size_t waited = 0;
	while (!queue.try_enqueue(token, write) && !queue.try_enqueue(token, write))
	{
		if (queue.size_approx() < 100 && queue.enqueue(token, write)) break;
		std::this_thread::sleep_for(std::chrono::milliseconds(10));
		waited++;
		if (waited >= 10)
		{
			if (queue.enqueue(token, write)) break;
		}
	}
}

void writeGAMToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, moodycamel::ConcurrentQueue<std::string*>& alignmentsOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	::google::protobuf::io::ZeroCopyOutputStream *raw_out = new ::google::protobuf::io::OstreamOutputStream(&strstr);
	::google::protobuf::io::GzipOutputStream *gzip_out = new ::google::protobuf::io::GzipOutputStream(raw_out);
	::google::protobuf::io::CodedOutputStream *coded_out = new ::google::protobuf::io::CodedOutputStream(gzip_out);
	coded_out->WriteVarint64(alignments.alignments.size());
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].alignment != nullptr);
		std::string s;
		alignments.alignments[i].alignment->SerializeToString(&s);
		coded_out->WriteVarint32(s.size());
		coded_out->WriteRaw(s.data(), s.size());
	}
	delete coded_out;
	delete gzip_out;
	delete raw_out;
	QueueInsertSlowly(token, alignmentsOut, strstr.str());
}

void writeJSONToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, moodycamel::ConcurrentQueue<std::string*>& alignmentsOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	google::protobuf::util::JsonPrintOptions options;
	options.preserve_proto_field_names = true;
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].alignment != nullptr);
		std::string s;
		google::protobuf::util::MessageToJsonString(*alignments.alignments[i].alignment, &s, options);
		strstr << s;
		strstr << '\n';
	}
	QueueInsertSlowly(token, alignmentsOut, strstr.str());
}

void writeGAFToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, moodycamel::ConcurrentQueue<std::string*>& alignmentsOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].GAFline.size() > 0);
		strstr << alignments.alignments[i].GAFline;
		strstr << '\n';
	}
	QueueInsertSlowly(token, alignmentsOut, strstr.str());
}

void writeCorrectedToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, const std::string& readName, const std::string& original, size_t maxOverlap, moodycamel::ConcurrentQueue<std::string*>& correctedOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	zstr::ostream *compressed;
	if (params.compressCorrected)
	{
		compressed = new zstr::ostream(strstr);
	}
	std::vector<Correction> corrections;
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].corrected.size() > 0);
		corrections.emplace_back();
		corrections.back().startIndex = alignments.alignments[i].alignmentStart;
		corrections.back().endIndex = alignments.alignments[i].alignmentEnd;
		corrections.back().corrected = alignments.alignments[i].corrected;
	}
	std::string corrected = getCorrected(original, corrections, maxOverlap);
	if (params.compressCorrected)
	{
		(*compressed) << ">" << readName << std::endl;
		(*compressed) << corrected << std::endl;
		delete compressed;
	}
	else
	{
		strstr << ">" << readName << std::endl;
		strstr << corrected << std::endl;
	}
	QueueInsertSlowly(token, correctedOut, strstr.str());
}

void writeCorrectedClippedToQueue(moodycamel::ProducerToken& token, const AlignerParams& params, moodycamel::ConcurrentQueue<std::string*>& correctedClippedOut, const AlignmentResult& alignments)
{
	std::stringstream strstr;
	zstr::ostream *compressed;
	if (params.compressClipped)
	{
		compressed = new zstr::ostream(strstr);
	}
	for (size_t i = 0; i < alignments.alignments.size(); i++)
	{
		assert(!alignments.alignments[i].alignmentFailed());
		assert(alignments.alignments[i].corrected.size() > 0);
		if (params.compressClipped)
		{
			(*compressed) << ">" << alignments.readName << "_" << i << "_" << alignments.alignments[i].alignmentStart << "_" << alignments.alignments[i].alignmentEnd << std::endl;
			(*compressed) << alignments.alignments[i].corrected << std::endl;
		}
		else
		{
			strstr << ">" << alignments.readName << "_" << i << "_" << alignments.alignments[i].alignmentStart << "_" << alignments.alignments[i].alignmentEnd << std::endl;
			strstr << alignments.alignments[i].corrected << std::endl;
		}
	}
	if (params.compressClipped)
	{
		delete compressed;
	}
	QueueInsertSlowly(token, correctedClippedOut, strstr.str());
}

std::vector<AlignmentGraph::MatrixPosition> traceToPoses(const AlignmentGraph& alignmentGraph, AlignmentResult::AlignmentItem &aln) {
	std::vector<AlignmentGraph::MatrixPosition> ret;
	auto trace = aln.trace->trace;
	size_t lastNode, lastOffset, lastLength;
	for (size_t j = 0; j < trace.size(); j++) {
		AlignmentGraph::MatrixPosition p = trace[j].DPposition;
		p.node = alignmentGraph.GetUnitigNode(p.node, p.nodeOffset);
		p.nodeOffset -= alignmentGraph.NodeOffset(p.node);
		if (j == 0) {
			lastNode = p.node;
			lastOffset = p.nodeOffset;
			lastLength = alignmentGraph.NodeLength(p.node);
			ret.emplace_back(lastNode, lastOffset, 0);
			lastOffset++;
		}
		else {
			if (p.node != lastNode) {
				while (lastOffset < lastLength) {
					ret.emplace_back(lastNode, lastOffset, 0);
					lastOffset++;
				}
				lastNode = p.node;
				lastLength = alignmentGraph.NodeLength(p.node);
				lastOffset = 0;
			}
			while (lastOffset <= p.nodeOffset) {
				ret.emplace_back(lastNode, lastOffset, 0);
				lastOffset++;
			}
		}
	}
	return ret;
}
std::vector<AlignmentGraph::MatrixPosition> pathToTrace(const AlignmentGraph& alignmentGraph, const std::vector<size_t> &path, size_t firstNodeOffset, size_t lastNodeOffset) {
	std::vector<AlignmentGraph::MatrixPosition> ret;
	for (size_t node : path) {
		size_t S = 0, L = alignmentGraph.NodeLength(node);
		if (node == path[0])
			S = firstNodeOffset;
		else if (node == path.back())
			L = lastNodeOffset + 1;
		AlignmentGraph::MatrixPosition p(node, S, 0);
		while (p.nodeOffset < L) {
			ret.push_back(p);
			p.nodeOffset++;
		}
	}
	return ret;
}
std::string traceToSequence(const AlignmentGraph& alignmentGraph, AlignmentResult::AlignmentItem &aln) {
	std::string ret = "";
	for (const AlignmentGraph::MatrixPosition &p : traceToPoses(alignmentGraph, aln)) 
		ret.push_back(alignmentGraph.NodeSequences(p.node, p.nodeOffset));
	// auto trace = aln.trace->trace;
	// size_t lastNode, lastOffset, lastLength;
	// for (size_t j = 0; j < trace.size(); j++) {
	// 	AlignmentGraph::MatrixPosition p = trace[j].DPposition;
	// 	p.node = alignmentGraph.GetUnitigNode(p.node, p.nodeOffset);
	// 	p.nodeOffset -= alignmentGraph.NodeOffset(p.node);
	// 	if (j == 0) {
	// 		lastNode = p.node;
	// 		lastOffset = p.nodeOffset;
	// 		lastLength = alignmentGraph.NodeLength(p.node);
	// 		ret.push_back(alignmentGraph.NodeSequences(lastNode, lastOffset));
	// 		lastOffset++;
	// 	}
	// 	else {
	// 		if (p.node != lastNode) {
	// 			while (lastOffset < lastLength) {
	// 				ret.push_back(alignmentGraph.NodeSequences(lastNode, lastOffset));
	// 				lastOffset++;
	// 			}
	// 			lastNode = p.node;
	// 			lastLength = alignmentGraph.NodeLength(p.node);
	// 			lastOffset = 0;
	// 		}
	// 		while (lastOffset <= p.nodeOffset) {
	// 			ret.push_back(alignmentGraph.NodeSequences(lastNode, lastOffset));
	// 			lastOffset++;
	// 		}
	// 	}
	// }
	return ret;
}

void runComponentMappings(const AlignmentGraph& alignmentGraph, moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>>& readFastqsQueue, std::atomic<bool>& readStreamingFinished, int threadnum, const Seeder& seeder, AlignerParams params, moodycamel::ConcurrentQueue<std::string*>& GAMOut, moodycamel::ConcurrentQueue<std::string*>& JSONOut, moodycamel::ConcurrentQueue<std::string*>& GAFOut, moodycamel::ConcurrentQueue<std::string*>& correctedOut, moodycamel::ConcurrentQueue<std::string*>& correctedClippedOut, moodycamel::ConcurrentQueue<std::string*>& deallocqueue, AlignmentStats& stats)
{
	moodycamel::ProducerToken GAMToken { GAMOut };
	moodycamel::ProducerToken JSONToken { JSONOut };
	moodycamel::ProducerToken GAFToken { GAFOut };
	moodycamel::ProducerToken correctedToken { correctedOut };
	moodycamel::ProducerToken clippedToken { correctedClippedOut };
	assertSetNoRead("Before any read");
	GraphAlignerCommon<size_t, int32_t, uint64_t>::AlignerGraphsizedState reusableState { alignmentGraph, std::max(params.initialBandwidth, params.rampBandwidth), !params.highMemory };
	AlignmentSelection::SelectionOptions selectionOptions;
	selectionOptions.method = params.alignmentSelectionMethod;
	selectionOptions.graphSize = alignmentGraph.SizeInBP();
	selectionOptions.ECutoff = params.selectionECutoff;
	if (params.preciseClipping)
	{
		selectionOptions.EValueCalc = EValueCalculator { params.preciseClippingIdentityCutoff };
	}
	else
	{
		// default 70% min identity threshold
		selectionOptions.EValueCalc = EValueCalculator { .7 };
	}
	BufferedWriter cerroutput;
	BufferedWriter coutoutput;
	if (params.verboseMode)
	{
		cerroutput = {std::cerr};
		coutoutput = {std::cout};
	}
	else if (params.shortVerboseMode)
		cerroutput = {std::cerr};
	while (true)
	{
		std::string* dealloc;
		while (deallocqueue.try_dequeue(dealloc))
		{
			delete dealloc;
		}
		std::shared_ptr<FastQ> fastq = nullptr;
		while (fastq == nullptr && !readFastqsQueue.try_dequeue(fastq))
		{
			bool tryBreaking = readStreamingFinished;
			if (!readFastqsQueue.try_dequeue(fastq) && tryBreaking) break;
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}
		if (fastq == nullptr) break;
		assertSetNoRead(fastq->seq_id);
		// cerroutput << "eds " << fastq->seq_id << " size " << fastq->sequence.size() << "bp" <<BufferedWriter::Flush;;
		selectionOptions.readSize = fastq->sequence.size();
		stats.reads += 1;
		stats.bpInReads += fastq->sequence.size();
		// auto AlllStart = std::chrono::system_clock::now();

		size_t tmpidx = stats.reads;
		std::string short_id;
		for (char c : fastq->seq_id)
			if (isspace(c))
				break;
			else
				short_id += c;
		// static std::string target_id = "S1_299";
		// if (!target_id.empty() && short_id != target_id) continue;

		// cerroutput << tmp << "  " << fastq->seq_id << " : " << fastq->sequence.length() << BufferedWriter::Flush;
		AlignmentResult alignments;

		size_t alntimems = 0;
		size_t clustertimems = 0;
		bool cont = false;
		// auto align_fn = [&seeder, &alignments, &reusableState, &alntimems, &clustertimems, &stats, &coutoutput, &cerroutput, &error, &params, &alignmentGraph] (const std::string &seq_id, const std::string &sequence) {
		auto align_fn = [&] (const std::string &seq_id, const std::string &sequence) {
			AlignmentResult alignments;
			try
			{
				if (seeder.mode != Seeder::Mode::None)
				{
					auto timeStart = std::chrono::system_clock::now();
					std::vector<SeedHit> seeds = seeder.getSeeds(seq_id, sequence);
					auto timeEnd = std::chrono::system_clock::now();
					size_t time = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart).count();
					coutoutput << "Read " << seq_id << " seeding took " << time << "ms" << BufferedWriter::Flush;
					stats.seeds += seeds.size();
					if (seeds.size() == 0)
					{
						coutoutput << "Read " << seq_id << " has no seed hits" << BufferedWriter::Flush;
						cerroutput << "Read " << seq_id << " has no seed hits" << BufferedWriter::Flush;
						coutoutput << "Read " << seq_id << " alignment failed" << BufferedWriter::Flush;
						cerroutput << "Read " << seq_id << " alignment failed" << BufferedWriter::Flush;
						if (params.outputCorrectedFile != "") writeCorrectedToQueue(correctedToken, params, seq_id, sequence, alignmentGraph.getDBGoverlap(), correctedOut, alignments);
						cont = true;
					}
					else {
						stats.seedsFound += seeds.size();
						stats.readsWithASeed += 1;
						stats.bpInReadsWithASeed += sequence.size();
						auto clusterTimeStart = std::chrono::system_clock::now();
						// if (params.colinearChaining)
						// 	OrderSeedsCLC(alignmentGraph, seeds);
						// else
							OrderSeeds(alignmentGraph, seeds);
						auto clusterTimeEnd = std::chrono::system_clock::now();
						clustertimems = std::chrono::duration_cast<std::chrono::milliseconds>(clusterTimeEnd - clusterTimeStart).count();
						coutoutput << "Read " << seq_id << " clustering took " << clustertimems << "ms" << BufferedWriter::Flush;
						auto alntimeStart = std::chrono::system_clock::now();
						alignments = AlignOneWay(alignmentGraph, seq_id, sequence, params.initialBandwidth, params.rampBandwidth, params.maxCellsPerSlice, !params.verboseMode, !params.tryAllSeeds, seeds, reusableState, !params.highMemory, params.forceGlobal, params.preciseClipping, params.seedClusterMinSize, params.seedExtendDensity, params.nondeterministicOptimizations, params.preciseClippingIdentityCutoff, params.Xdropcutoff);
						auto alntimeEnd = std::chrono::system_clock::now();
						alntimems = std::chrono::duration_cast<std::chrono::milliseconds>(alntimeEnd - alntimeStart).count();
					}
				}
				else if (params.optimalDijkstra)
				{
					auto alntimeStart = std::chrono::system_clock::now();
					alignments = AlignOneWayDijkstra(alignmentGraph, seq_id, sequence, !params.verboseMode, reusableState, params.forceGlobal, params.preciseClipping);
					auto alntimeEnd = std::chrono::system_clock::now();
					alntimems = std::chrono::duration_cast<std::chrono::milliseconds>(alntimeEnd - alntimeStart).count();
				}
				else
				{
					auto alntimeStart = std::chrono::system_clock::now();
					alignments = AlignOneWay(alignmentGraph, seq_id, sequence, params.initialBandwidth, params.rampBandwidth, !params.verboseMode, reusableState, !params.highMemory, params.forceGlobal, params.preciseClipping, params.nondeterministicOptimizations, params.preciseClippingIdentityCutoff, params.Xdropcutoff, params.DPRestartStride);
					auto alntimeEnd = std::chrono::system_clock::now();
					alntimems = std::chrono::duration_cast<std::chrono::milliseconds>(alntimeEnd - alntimeStart).count();
				}
			}
			catch (const ThreadReadAssertion::AssertionFailure& a)
			{
				coutoutput << "Read " << seq_id << " alignment failed (assertion!)" << BufferedWriter::Flush;
				cerroutput << "Read " << seq_id << " alignment failed (assertion!)" << BufferedWriter::Flush;
				reusableState.clear();
				stats.assertionBroke = true;
				cont = true;
			}
			return alignments;
		};
		
		if (!params.colinearChaining) {
			alignments = align_fn(fastq->seq_id, fastq->sequence);
			if (cont)
				continue;
		}
		if (params.colinearChaining) {
			// check whether colinear is necessary
			bool necessary = true;
			AlignmentResult long_alignments;
			size_t long_edit_distance;
			// long_alignments = align_fn(fastq->seq_id, fastq->sequence);
			// if (long_alignments.alignments.empty())
			// 	necessary = true;
			// else {
			// 	std::string long_pathseq = traceToSequence(alignmentGraph, long_alignments.alignments[0]);
			// 	if (long_pathseq.length() < 0.9 * fastq->sequence.length())
			// 		necessary = true;
			// 	else {
			// 		EdlibAlignResult result = edlibAlign(long_pathseq.c_str(), long_pathseq.length(), fastq->sequence.c_str(), fastq->sequence.length(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
			// 		if (result.status != EDLIB_STATUS_OK) {
			// 			necessary = true;
			// 			long_edit_distance = fastq->sequence.length();
			// 		}
			// 		else {
			// 			long_edit_distance = result.editDistance;
			// 			if (long_edit_distance > 0.2 * fastq->sequence.length())
			// 				necessary = true;
			// 		}
			// 		edlibFreeAlignResult(result);
			// 	}
			// }
			if (!necessary)
				alignments = std::move(long_alignments);
			else if (necessary) {
				// start spliting long reads into short reads for anchors
				std::vector<AlignmentGraph::Anchor> A;
				std::vector<std::vector<GraphAlignerCommon<size_t, int32_t, uint64_t>::TraceItem>> Apos;
				std::vector<SeedHit> seeds = seeder.getSeeds(fastq->seq_id, fastq->sequence);
				if (seeds.size() == 0)
					continue;
				stats.seedsFound += seeds.size();
				stats.readsWithASeed += 1;
				stats.bpInReadsWithASeed += fastq->sequence.size();
				OrderSeeds(alignmentGraph, seeds);
				std::sort(seeds.begin(), seeds.end(), [](const SeedHit& left, const SeedHit& right) { return left.seqPos < right.seqPos; });
				size_t len = params.colinearSplitLen, sep = params.colinearSplitGap; // by default 150,33 
				size_t sl = 0, sr = 0;
				
				auto anchorsStart = std::chrono::system_clock::now();
				for (size_t l = 0; l + len <= fastq->sequence.length(); l += sep) {
					while (sr < seeds.size() && seeds[sr].seqPos + seeds[sr].matchLen <= l + len)
						sr++;
					while (sl < sr && seeds[sl].seqPos < l)
						sl++;
					// cerroutput << short_id << " : " << l << " / " << fastq->sequence.length() << " " << sl << " " << sr << BufferedWriter::Flush;
					if (sl >= sr)
						continue;
					std::string seq = fastq->sequence.substr(l, len);
					std::string name = short_id + "_" + std::to_string(l) + "_" + std::to_string(l + len - 1);
					try {
						// for (size_t k = sl; k < sr; k++)
						// 	cerroutput << k << " : " << seeds[k].seqPos << BufferedWriter::Flush;
						// std::vector<SeedHit> tmp_seeds;
						// for (size_t k = sl; k < sr; k++) {
						// 	tmp_seeds.push_back(seeds[k]);
						// 	tmp_seeds.back().seqPos -= l;
						// }
						if (seeder.mode != Seeder::Mode::None) {
							alignments = AlignOneWay(alignmentGraph, name, seq, params.initialBandwidth, params.rampBandwidth, params.maxCellsPerSlice, !params.verboseMode, !params.tryAllSeeds, seeds, reusableState, !params.highMemory, params.forceGlobal, params.preciseClipping, params.seedClusterMinSize, params.seedExtendDensity, params.nondeterministicOptimizations, params.preciseClippingIdentityCutoff, params.Xdropcutoff, sl, sr, l);
							// alignments = AlignOneWay(alignmentGraph, name, seq, params.initialBandwidth, params.rampBandwidth, params.maxCellsPerSlice, !params.verboseMode, !params.tryAllSeeds, tmp_seeds, reusableState, !params.highMemory, params.forceGlobal, params.preciseClipping, params.seedClusterMinSize, params.seedExtendDensity, params.nondeterministicOptimizations, params.preciseClippingIdentityCutoff, params.Xdropcutoff, 0, tmp_seeds.size(), 0);
						}
					}
					catch (const ThreadReadAssertion::AssertionFailure& a) {
						coutoutput << "Read " << short_id << " alignment failed (assertion!)" << BufferedWriter::Flush;
						cerroutput << "Read " << short_id << " alignment failed (assertion!)" << BufferedWriter::Flush;
						reusableState.clear();
						stats.assertionBroke = true;
						cont = true;
					}
					if (cont)
						continue;
					// all short alignments are used as anchors
					stats.seedsExtended += alignments.seedsExtended;
					for (size_t i = 0; i < alignments.alignments.size(); i++) {
						AlignmentGraph::Anchor anchor = {{}, l, l + len - 1};
						AlignmentResult::AlignmentItem& alignment = alignments.alignments[i];
						if (alignment.alignmentFailed())
							continue;
						auto trace = alignment.trace->trace;
						if (trace.size() == 0)
							continue;
						for (size_t j = 0; j < trace.size(); j++) {
							size_t node = trace[j].DPposition.node;
							size_t nodeOffset = trace[j].DPposition.nodeOffset;
							node = alignmentGraph.GetUnitigNode(node, nodeOffset);
							if (anchor.path.empty() || node != anchor.path.back())
								anchor.path.push_back(node);
						}
						A.push_back(anchor);
						Apos.push_back({ trace[0], trace.back() });
						for (size_t j = 0; j < Apos.back().size(); j++) {
							AlignmentGraph::MatrixPosition &p = Apos.back()[j].DPposition;
							p.seqPos += l;
							p.node = alignmentGraph.GetUnitigNode(p.node, p.nodeOffset);
							p.nodeOffset -= alignmentGraph.NodeOffset(p.node);
						}
					}
				}
				auto anchorsEnd = std::chrono::system_clock::now();
				auto anchorsms = std::chrono::duration_cast<std::chrono::milliseconds>(anchorsEnd - anchorsStart).count();
				// cerroutput << short_id << " : chaining " << A.size() << " anchors" << BufferedWriter::Flush;
				auto clcStart = std::chrono::system_clock::now();
				std::vector<size_t> ids = alignmentGraph.colinearChaining(A, params.colinearGap);
				auto clcEnd = std::chrono::system_clock::now();
				auto clcms = std::chrono::duration_cast<std::chrono::milliseconds>(clcEnd - clcStart).count();
				alignments.alignments.clear();
				GraphAlignerCommon<size_t, int32_t, uint64_t>::OnewayTrace trace;
				std::vector<AlignmentGraph::MatrixPosition> longest, tmp;
				std::vector<size_t> pos_path;
				std::unordered_set<size_t> nodes;
				size_t firstNodeOffset, lastNodeOffset;
				auto connectStart = std::chrono::system_clock::now();
				int ai_idx = 0;
				// cerroutput << short_id << " : chaining " << ids.size() << " / " << A.size() << " anchors" << BufferedWriter::Flush;
				int one_node_overlaps_tmp = 0, one_node_overlaps_now = 0, one_node_overlaps_all = 0;
				for (size_t ai : ids) {
					for (size_t j =0 ; j<Apos[ai].size();j++) {
						AlignmentGraph::MatrixPosition &p = Apos[ai][j].DPposition;
							// if (alignmentGraph.NodeLength(p.node) <= p.nodeOffset) std::cerr << "????" << ai << " :" << j << std::endl;
					}
				}
				for (size_t ai : ids) {
					const AlignmentGraph::Anchor &anchor = A[ai];
					// std::cerr << short_id << " : connect " << ai << " " << (ai_idx++) << "/" << ids.size() << "  " << pos_path.size() << "   longsest" << longest.size()<< std::endl;
					// if (ai_idx>1) std::cerr << "ccc " << alignmentGraph.getChainPath(A[ids[ai_idx-2]].path.back(), A[ids[ai_idx-1]].path[0], -1).size() << std::endl;
					// std::cerr << "now " << Apos[ai][0].DPposition.nodeOffset << "  " << Apos[ai].back().DPposition.nodeOffset << std::endl;
					// std::cerr << "back " << Apos[ai].back().DPposition.node << " " << alignmentGraph.NodeLength(Apos[ai].back().DPposition.node) << std::endl;
					if (pos_path.empty()) {
						pos_path = anchor.path;
						firstNodeOffset = Apos[ai][0].DPposition.nodeOffset;
						lastNodeOffset = Apos[ai].back().DPposition.nodeOffset;
						for (size_t j : pos_path)
							nodes.insert(j);
					}
					else {
						if (anchor.path[0] == pos_path.back()) {
							one_node_overlaps_tmp++;
							one_node_overlaps_all++;
						}
						bool gap = anchor.path[0] == pos_path.back() && params.colinearGap != -1 && (long long)Apos[ai][0].DPposition.nodeOffset - (long long)lastNodeOffset > params.colinearGap + 1;
						if (gap) {
							// std::cerr << " ? " << anchor.path[0] << " to " << pos_path.back() << " : " << Apos[ai][0].DPposition.nodeOffset << " " << lastNodeOffset << " || " << Apos[ai][0].DPposition.node << std::endl;
						}
						std::vector<size_t> path;
						if (!nodes.count(anchor.path[0]) && pos_path.back() != Apos[ai][0].DPposition.node) {
							long long gapLimit = params.colinearGap;
							if (gapLimit != -1)
								gapLimit -= (long long)Apos[ai][0].DPposition.nodeOffset + (long long)(alignmentGraph.NodeLength(pos_path.back()) - (long long)lastNodeOffset - 1);
							path = alignmentGraph.getChainPath(pos_path.back(), Apos[ai][0].DPposition.node, gapLimit);
							// std::cerr << "gap _bfs_ " << (int)(gap) << " " << path.size() << "  " << gapLimit << std::endl;
							if (path.empty())
								gap = true;
						}
						if (gap) {
							// std::cerr << "gap _gap_ " << (int)(gap) << " " << path.size() << "  " << std::endl;
							tmp = pathToTrace(alignmentGraph, pos_path, firstNodeOffset, lastNodeOffset);
				// cerroutput << tmp.size() << " : " << pos_path.size() << " at " << firstNodeOffset << " " << lastNodeOffset << BufferedWriter::Flush;
							if (longest.size() < tmp.size()) {
								longest.swap(tmp);
								one_node_overlaps_now = one_node_overlaps_tmp;
							}
							nodes.clear();
							pos_path.clear();
							firstNodeOffset = Apos[ai][0].DPposition.nodeOffset;
						}
						else
							for (size_t j : path) 
								if (!nodes.count(j)) {
									nodes.insert(j);
									pos_path.push_back(j);
								}
						// for (size_t j : anchor.path) std::cerr << " " << j << (int)(nodes.count(j));
						//  std::cerr << std::endl;
						for (size_t j : anchor.path) 
							if (!nodes.count(j)) {
								nodes.insert(j);
								pos_path.push_back(j);
							}
						lastNodeOffset = Apos[ai].back().DPposition.nodeOffset;
					}
				}
				if (!pos_path.empty()) {
					tmp = pathToTrace(alignmentGraph, pos_path, firstNodeOffset, lastNodeOffset);
					// cerroutput << tmp.size() << " : " << pos_path.size() << " at " << firstNodeOffset << " " << lastNodeOffset << BufferedWriter::Flush;
					// cerroutput << "last " << pos_path.back() <<" " << lastNodeOffset << alignmentGraph.NodeLength(pos_path.back()) << BufferedWriter::Flush;
					if (longest.size() < tmp.size()) {
						longest.swap(tmp);
						one_node_overlaps_now = one_node_overlaps_tmp;
					}
				}
				// cerroutput << short_id << " : chaining " << longest.size() << " bps" << BufferedWriter::Flush;

				std::string pathseq = "";
				// not convert back to original node ids, to call alignmentGraph.NodeSequences()
				for (AlignmentGraph::MatrixPosition &p : longest) {
				// cerroutput << p.node << " : " << p.nodeOffset << " bps " << alignmentGraph.NodeLength(p.node) << BufferedWriter::Flush;

					pathseq.push_back(alignmentGraph.NodeSequences(p.node, p.nodeOffset));
				}
				// align by edit distance to get trace, using edlib
				size_t alnScore = 0;
				if (params.fastMode) {
					if (ids.size() > 0) {
						size_t x = A[ids[0]].x, y = A[ids.back()].y;
						for (size_t j = 0; j < longest.size(); j++) {
							longest[j].seqPos = std::min(y, x + j);
							if (alignmentGraph.NodeSequences(longest[j].node, longest[j].nodeOffset) != fastq->sequence[longest[j].seqPos])
								alnScore++;
						}
					}
				}
				else {
					EdlibAlignResult result = edlibAlign(pathseq.c_str(), pathseq.length(), fastq->sequence.c_str(), fastq->sequence.length(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
					if (result.status != EDLIB_STATUS_OK)
						longest.clear();
					else {
						alnScore = result.editDistance;
						// std::cerr << short_id << ": eds " << result.editDistance << " from  len=" << pathseq.length()<<" " << fastq->sequence.length()<< std::endl;
						std::vector<AlignmentGraph::MatrixPosition> trace;
						trace.reserve(result.alignmentLength);
						size_t start = result.startLocations[0], end = result.endLocations[0];
						size_t pos_i = 0, seq_i = result.startLocations[0];
						for (size_t j = 0; j < result.alignmentLength; j++) {
							AlignmentGraph::MatrixPosition p(longest[pos_i].node, longest[pos_i].nodeOffset, seq_i);
							// char a = alignmentGraph.NodeSequences(p.node, p.nodeOffset);
							// char b = fastq->sequence[seq_i];
							trace.push_back(p);
							unsigned char c = result.alignment[j];

							if (c == 0 || c == 3) { // match or mismatch
								pos_i++;
								seq_i++;
							}
							else if (c == 1)
								pos_i++;
							else if (c == 2)
								seq_i++;
							seq_i = std::min(seq_i, fastq->sequence.length() - 1);
							pos_i = std::min(pos_i, longest.size() - 1);
							
						}
						// cerroutput << "eds " << trace.size() << " " << longest.size() << " : " << pos_i << " / " << seq_i << BufferedWriter::Flush;;
						longest.swap(trace);
					}
					edlibFreeAlignResult(result);
				}
				
				for (size_t i = 0; i < longest.size(); i++) {
					bool nodeSwitch = false;
					if (i + 1 < longest.size() && longest[i].node != longest[i + 1].node)
						nodeSwitch = true;
					trace.trace.emplace_back(longest[i], nodeSwitch, fastq->sequence, alignmentGraph);
					AlignmentGraph::MatrixPosition &p = trace.trace.back().DPposition;
					p.nodeOffset += alignmentGraph.NodeOffset(p.node);
					p.node = alignmentGraph.NodeID(p.node);
				}
				
				if (trace.trace.size() > 0) {
					AlignmentResult::AlignmentItem result { std::move(trace), 0, std::numeric_limits<size_t>::max() };
					assert(result.trace->trace.size() > 0);
					result.alignmentScore = alnScore; //result.trace->score;
					result.alignmentStart = result.trace->trace[0].DPposition.seqPos;
					result.alignmentEnd = result.trace->trace.back().DPposition.seqPos + 1;
					alignments.alignments.push_back(result);
				}
				auto connectEnd = std::chrono::system_clock::now();
				auto connectms = std::chrono::duration_cast<std::chrono::milliseconds>(connectEnd - connectStart).count();
				// bool better = (long_alignments.alignments.empty() || long_edit_distance > alignments.alignments.back().alignmentScore);

				if (params.shortVerboseMode || params.verboseMode)
					cerroutput << tmpidx << " " << short_id << " len=" << fastq->sequence.length() << " : "
					<< "chained " << ids.size() << " / " << A.size() << " anchors, actual " << longest.size() << " bps, "
					<< "time " << anchorsms << " " << clcms << " " << connectms << "  "
					<< "score=" << alnScore// << " better? " << (better?"Yes":"No") 
					<< " one_node_overlaps=" << one_node_overlaps_now << " / " << one_node_overlaps_all
					<< BufferedWriter::Flush;

				// compare alignments
				// if (!better) {
				// 	alignments = std::move(long_alignments);
				// }
			}
		}

		stats.allAlignmentsCount += alignments.alignments.size();

		coutoutput << "Read " << fastq->seq_id << " alignment took " << alntimems << "ms" << BufferedWriter::Flush;
		if (alignments.alignments.size() > 0) alignments.alignments = AlignmentSelection::SelectAlignments(alignments.alignments, selectionOptions);

		// for (size_t i = 0; i < alignments.alignments.size(); i++) {
		// 	// c++ helloWorld.cpp edlib/src/edlib.cpp -o helloWorld -I edlib/include.
		// 	std::string pathseq = traceToSequence(alignmentGraph, alignments.alignments[i]);
		// 	// global edit distance
		// 	EdlibAlignResult result = edlibAlign(pathseq.c_str(), pathseq.length(), fastq->sequence.c_str(), fastq->sequence.length(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		// 	// local edit distance
		// 	EdlibAlignResult result2 = edlibAlign(pathseq.c_str(), pathseq.length(), fastq->sequence.c_str(), fastq->sequence.length(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
		// 	// if (result.status == EDLIB_STATUS_OK) {
		// 	// 	printf("edit_distance('hello', 'world!') = %d\n", result.editDistance);
		// 	// }
		// 	// std::cerr << i << "   " << fastq->sequence.length()<< " : " << alignments.alignments[i].alignmentScore << "   ---   " << result.editDistance << "   ---   " << result2.editDistance << std::endl;

		// 	std::vector<AlignmentGraph::MatrixPosition> poss = traceToPoses(alignmentGraph, alignments.alignments[i]), trace;
		// 	trace.reserve(result.alignmentLength);
		// 	size_t start = result.startLocations[0], end = result.endLocations[0];
		// 	size_t pos_i = 0, seq_i = result.startLocations[0];
		// 	for (size_t j = 0; j < result.alignmentLength; j++) {
		// 		AlignmentGraph::MatrixPosition p(poss[pos_i].node, poss[pos_i].nodeOffset, seq_i);
		// 		char a = alignmentGraph.NodeSequences(p.node, p.nodeOffset);
		// 		char b = fastq->sequence[seq_i];
		// 		p.nodeOffset += alignmentGraph.NodeOffset(p.node);
		// 		p.node = alignmentGraph.NodeID(p.node);
		// 		trace.push_back(p);
		// 		unsigned char c = result.alignment[j];
		// 		if (c == 0 || c == 3) { // match or mismatch
		// 			pos_i++;
		// 			seq_i++;
		// 			// if (c == 0 && a != b) std::cerr << a << " " << b << " " << c << std::endl;
		// 			// if (c == 1 && a == b) std::cerr << a << " " << b << " " << c << std::endl;
		// 		}
		// 		else if (c == 1)
		// 			pos_i++;
		// 		else if (c == 2)
		// 			seq_i++;
		// 	}
		// 	// for (size_t j = 0; j < result.numLocations && j < 1; j++) {
		// 	// 	std::cerr << " ? " << j << " " << result.startLocations[j] << " : " << result.endLocations[j] << std::endl;
		// 	// }
		// 	// std::cerr << result.numLocations << " " << result.alignmentLength << " @ " << pos_i << " / " << poss.size() << " " << seq_i << " / " << fastq->sequence.length() << std::endl;

		// 	edlibFreeAlignResult(result);
		// 	edlibFreeAlignResult(result2);
		// }

		//failed alignment, don't output
		if (alignments.alignments.size() == 0)
		{
			coutoutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			cerroutput << "Read " << fastq->seq_id << " alignment failed" << BufferedWriter::Flush;
			try
			{
				if (params.outputCorrectedFile != "") writeCorrectedToQueue(correctedToken, params, fastq->seq_id, fastq->sequence, alignmentGraph.getDBGoverlap(), correctedOut, alignments);
			}
			catch (const ThreadReadAssertion::AssertionFailure& a)
			{
				reusableState.clear();
				stats.assertionBroke = true;
				continue;
			}
			continue;
		}

		stats.seedsExtended += alignments.seedsExtended;
		stats.readsWithAnAlignment += 1;

		size_t totalcells = 0;
		for (size_t i = 0; i < alignments.alignments.size(); i++)
		{
			totalcells += alignments.alignments[i].cellsProcessed;
		}
		
		std::sort(alignments.alignments.begin(), alignments.alignments.end(), [](const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right) { return left.alignmentStart < right.alignmentStart; });

		if (params.outputGAMFile != "" || params.outputJSONFile != "")
		{
			for (size_t i = 0; i < alignments.alignments.size(); i++)
			{
				AddAlignment(fastq->seq_id, fastq->sequence, alignments.alignments[i]);
				replaceDigraphNodeIdsWithOriginalNodeIds(*alignments.alignments[i].alignment, alignmentGraph);
			}
		}

		if (params.outputGAFFile != "")
		{
			for (size_t i = 0; i < alignments.alignments.size(); i++)
			{
				AddGAFLine(alignmentGraph, fastq->seq_id, fastq->sequence, alignments.alignments[i], params.cigarMatchMismatchMerge);
			}
		}
		
		std::sort(alignments.alignments.begin(), alignments.alignments.end(), [](const AlignmentResult::AlignmentItem& left, const AlignmentResult::AlignmentItem& right) { return left.alignmentStart < right.alignmentStart; });

		std::string alignmentpositions;

		for (size_t i = 0; i < alignments.alignments.size(); i++)
		{
			stats.alignments += 1;
			size_t alignmentSize = alignments.alignments[i].alignmentEnd - alignments.alignments[i].alignmentStart;
			if (alignmentSize == fastq->sequence.size())
			{
				stats.fullLengthAlignments += 1;
				stats.bpInFullAlignments += alignmentSize;
			}
			stats.bpInAlignments += alignmentSize;
			if (params.outputCorrectedFile != "" || params.outputCorrectedClippedFile != "") AddCorrected(alignments.alignments[i]);
			alignmentpositions += std::to_string(alignments.alignments[i].alignmentStart) + "-" + std::to_string(alignments.alignments[i].alignmentEnd) + ", ";
		}

		alignmentpositions.pop_back();
		alignmentpositions.pop_back();
		coutoutput << "Read " << fastq->seq_id << " aligned by thread " << threadnum << " with positions: " << alignmentpositions << " (read " << fastq->sequence.size() << "bp)" << BufferedWriter::Flush;

		try
		{
			if (params.outputGAMFile != "") writeGAMToQueue(GAMToken, params, GAMOut, alignments);
			if (params.outputJSONFile != "") writeJSONToQueue(JSONToken, params, JSONOut, alignments);
			if (params.outputGAFFile != "") writeGAFToQueue(GAFToken, params, GAFOut, alignments);
			if (params.outputCorrectedFile != "") writeCorrectedToQueue(correctedToken, params, fastq->seq_id, fastq->sequence, alignmentGraph.getDBGoverlap(), correctedOut, alignments);
			if (params.outputCorrectedClippedFile != "") writeCorrectedClippedToQueue(clippedToken, params, correctedClippedOut, alignments);
		}
		catch (const ThreadReadAssertion::AssertionFailure& a)
		{
			reusableState.clear();
			stats.assertionBroke = true;
			continue;
		}
		// auto Alllms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - AlllStart).count();
		// cerroutput << tmpidx << " " << short_id << " time " << Alllms << BufferedWriter::Flush;

	}
	assertSetNoRead("After all reads");
	coutoutput << "Thread " << threadnum << " finished" << BufferedWriter::Flush;
}

AlignmentGraph getGraph(std::string graphFile, MummerSeeder** mxmSeeder, const AlignerParams& params)
{
	bool loadMxmSeeder = params.mumCount > 0 || params.memCount > 0;
	if (is_file_exist(graphFile)){
		std::cout << "Load graph from " << graphFile << std::endl;
	}
	else{
		std::cerr << "No graph file exists" << std::endl;
		std::exit(0);
	}
	try
	{
		if (graphFile.substr(graphFile.size()-3) == ".vg")
		{
			if (loadMxmSeeder)
			{
				auto graph = CommonUtils::LoadVGGraph(graphFile);
				if (loadMxmSeeder)
				{
					std::cout << "Build MUM/MEM seeder from the graph" << std::endl;
					*mxmSeeder = new MummerSeeder { graph, params.seederCachePrefix };
				}
				std::cout << "Build alignment graph" << std::endl;
				auto result = DirectedGraph::BuildFromVG(graph);
				return result;
			}
			else
			{
				return DirectedGraph::StreamVGGraphFromFile(graphFile);
			}
		}
		else if (graphFile.substr(graphFile.size() - 4) == ".gfa")
		{
			auto graph = GfaGraph::LoadFromFile(graphFile, true);
			if (loadMxmSeeder)
			{
				std::cout << "Build MUM/MEM seeder from the graph" << std::endl;
				*mxmSeeder = new MummerSeeder { graph, params.seederCachePrefix };
			}
			std::cout << "Build alignment graph" << std::endl;
			auto result = DirectedGraph::BuildFromGFA(graph);
			return result;
		}
		else
		{
			std::cerr << "Unknown graph type (" << graphFile << ")" << std::endl;
			std::exit(0);
		}
	}
	catch (const CommonUtils::InvalidGraphException& e)
	{
		std::cout << "Error in the graph: " << e.what() << std::endl;
		std::cerr << "Error in the graph: " << e.what() << std::endl;
		std::exit(1);
	}
}

void alignReads(AlignerParams params)
{
	assertSetNoRead("Preprocessing");
	std::cout << "Co-linear chaining " << (params.colinearChaining ? "on" : "off");
	if (params.colinearChaining) {
		size_t len = params.colinearSplitLen, sep = params.colinearSplitGap, gap = params.colinearGap; // by default 150,33,1000
		std::cout << " splits=(" << len << "," << sep << "," << gap << ")";
	}
	std::cout << std::endl;

	const std::unordered_map<std::string, std::vector<SeedHit>>* seedHitsToThreads = nullptr;
	std::unordered_map<std::string, std::vector<SeedHit>> seedHits;
	MummerSeeder* mummerseeder = nullptr;
	auto alignmentGraph = getGraph(params.graphFile, &mummerseeder, params);
	if (params.generatePath) {
		std::vector<size_t> path = alignmentGraph.generatePath(params.fastqFiles[0], params.outputGAMFile, params.generatePathSeed);
		return;
	}
	if (params.graphStatistics) {
		alignmentGraph.buildMPC();
		return;
	}
	if (params.colinearChaining) {
		bool mpcReady = false;
		if (params.IndexMpcFile != "" && is_file_exist(params.IndexMpcFile)) {
			alignmentGraph.loadMPC(params.IndexMpcFile);
			mpcReady = true; //alignmentGraph.checkMPC();
		}
		if (!mpcReady) {
			alignmentGraph.buildMPC();
			alignmentGraph.saveMPC(params.IndexMpcFile);
		}
	}
	bool loadMinimizerSeeder = params.minimizerSeedDensity != 0;
	MinimizerSeeder* minimizerseeder = nullptr;
	if (loadMinimizerSeeder)
	{
		std::cout << "Build minimizer seeder from the graph" << std::endl;
		minimizerseeder = new MinimizerSeeder(alignmentGraph, params.minimizerLength, params.minimizerWindowSize, params.numThreads, 1.0 - params.minimizerDiscardMostNumerousFraction);
		if (!minimizerseeder->canSeed())
		{
			std::cout << "Warning: Minimizer seeder has no seed hits. Reads cannot be aligned. Try unchopping the graph with vg or a different seeding mode" << std::endl;
		}
	}

	if (params.seedFiles.size() > 0)
	{
		for (auto file : params.seedFiles)
		{
			if (is_file_exist(file)){
				std::cout << "Load seeds from " << file << std::endl;
				std::ifstream seedfile { file, std::ios::in | std::ios::binary };
				size_t numSeeds = 0;
				std::function<void(vg::Alignment&)> alignmentLambda = [&seedHits, &numSeeds](vg::Alignment& seedhit) {
					seedHits[seedhit.name()].emplace_back(seedhit.path().mapping(0).position().node_id(), seedhit.path().mapping(0).position().offset(), seedhit.query_position(), seedhit.path().mapping(0).edit(0).from_length(), seedhit.path().mapping(0).edit(0).from_length(), seedhit.path().mapping(0).position().is_reverse());
					numSeeds += 1;
				};
				stream::for_each(seedfile, alignmentLambda);
				std::cout << numSeeds << " seeds" << std::endl;
			}
			else {
				std::cerr << "No seeds file exists" << std::endl;
				std::exit(0);
			}
		}
		seedHitsToThreads = &seedHits;
	}

	Seeder seeder { params, seedHitsToThreads, mummerseeder, minimizerseeder };

	switch(seeder.mode)
	{
		case Seeder::Mode::File:
			std::cout << "Seeds from file" << std::endl;
			break;
		case Seeder::Mode::Mum:
			std::cout << "MUM seeds, min length " << seeder.mxmLength;
			if (seeder.mumCount != std::numeric_limits<size_t>::max()) std::cout << ", max count " << seeder.mumCount;
			std::cout << std::endl;
			break;
		case Seeder::Mode::Mem:
			std::cout << "MEM seeds, min length " << seeder.mxmLength;
			if (seeder.memCount != std::numeric_limits<size_t>::max()) std::cout << ", max count " << seeder.memCount;
			std::cout << std::endl;
			break;
		case Seeder::Mode::Minimizer:
			std::cout << "Minimizer seeds, length " << seeder.minimizerLength << ", window size " << seeder.minimizerWindowSize << ", density " << seeder.minimizerSeedDensity << std::endl;
			break;
		case Seeder::Mode::None:
			if (params.optimalDijkstra)
			{
				std::cout << "Optimal alignment. VERY SLOW!" << std::endl;
			}
			else
			{
				std::cout << "No seeds, calculate the entire first row. VERY SLOW!" << std::endl;
			}
			break;
	}
	if (seeder.mode != Seeder::Mode::None) std::cout << "Seed cluster size " << params.seedClusterMinSize << std::endl;
	if (seeder.mode != Seeder::Mode::None && params.seedExtendDensity != -1) std::cout << "Extend up to best " << params.seedExtendDensity << " fraction of seeds" << std::endl;

	if (!params.optimalDijkstra) std::cout << "Initial bandwidth " << params.initialBandwidth;
	if (params.rampBandwidth > 0) std::cout << ", ramp bandwidth " << params.rampBandwidth;
	if (params.maxCellsPerSlice != std::numeric_limits<size_t>::max()) std::cout << ", tangle effort " << params.maxCellsPerSlice;
	std::cout << std::endl;

	if (params.selectionECutoff != -1) std::cout << "Discard alignments with an E-value > " << params.selectionECutoff << std::endl;

	if (params.outputGAMFile != "") std::cout << "write alignments to " << params.outputGAMFile << std::endl;
	if (params.outputJSONFile != "") std::cout << "write alignments to " << params.outputJSONFile << std::endl;
	if (params.outputGAFFile != "") std::cout << "write alignments to " << params.outputGAFFile << std::endl;
	if (params.outputCorrectedFile != "") std::cout << "write corrected reads to " << params.outputCorrectedFile << std::endl;
	if (params.outputCorrectedClippedFile != "") std::cout << "write corrected & clipped reads to " << params.outputCorrectedClippedFile << std::endl;

	std::vector<std::thread> threads;

	assertSetNoRead("Running alignments");

	moodycamel::ConcurrentQueue<std::string*> outputGAM { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::string*> outputGAF { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::string*> outputJSON { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::string*> deallocAlns;
	moodycamel::ConcurrentQueue<std::string*> outputCorrected { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::string*> outputCorrectedClipped { 50, params.numThreads, params.numThreads };
	moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>> readFastqsQueue;
	std::atomic<bool> readStreamingFinished { false };
	std::atomic<bool> allThreadsDone { false };
	std::atomic<bool> GAMWriteDone { false };
	std::atomic<bool> GAFWriteDone { false };
	std::atomic<bool> JSONWriteDone { false };
	std::atomic<bool> correctedWriteDone { false };
	std::atomic<bool> correctedClippedWriteDone { false };

	std::cout << "Align" << std::endl;
	AlignmentStats stats;
	std::thread fastqThread { [files=params.fastqFiles, &readFastqsQueue, &readStreamingFinished]() { readFastqs(files, readFastqsQueue, readStreamingFinished); } };
	std::thread GAMwriterThread { [file=params.outputGAMFile, &outputGAM, &deallocAlns, &allThreadsDone, &GAMWriteDone, verboseMode=params.verboseMode]() { if (file != "") consumeBytesAndWrite(file, outputGAM, deallocAlns, allThreadsDone, GAMWriteDone, verboseMode, false); else GAMWriteDone = true; } };
	std::thread GAFwriterThread { [file=params.outputGAFFile, &outputGAF, &deallocAlns, &allThreadsDone, &GAFWriteDone, verboseMode=params.verboseMode]() { if (file != "") consumeBytesAndWrite(file, outputGAF, deallocAlns, allThreadsDone, GAFWriteDone, verboseMode, false); else GAFWriteDone = true; } };
	std::thread JSONwriterThread { [file=params.outputJSONFile, &outputJSON, &deallocAlns, &allThreadsDone, &JSONWriteDone, verboseMode=params.verboseMode]() { if (file != "") consumeBytesAndWrite(file, outputJSON, deallocAlns, allThreadsDone, JSONWriteDone, verboseMode, true); else JSONWriteDone = true; } };
	std::thread correctedWriterThread { [file=params.outputCorrectedFile, &outputCorrected, &deallocAlns, &allThreadsDone, &correctedWriteDone, verboseMode=params.verboseMode, uncompressed=!params.compressCorrected]() { if (file != "") consumeBytesAndWrite(file, outputCorrected, deallocAlns, allThreadsDone, correctedWriteDone, verboseMode, uncompressed); else correctedWriteDone = true; } };
	std::thread correctedClippedWriterThread { [file=params.outputCorrectedClippedFile, &outputCorrectedClipped, &deallocAlns, &allThreadsDone, &correctedClippedWriteDone, verboseMode=params.verboseMode, uncompressed=!params.compressClipped]() { if (file != "") consumeBytesAndWrite(file, outputCorrectedClipped, deallocAlns, allThreadsDone, correctedClippedWriteDone, verboseMode, uncompressed); else correctedClippedWriteDone = true; } };

	for (size_t i = 0; i < params.numThreads; i++)
	{
		threads.emplace_back([&alignmentGraph, &readFastqsQueue, &readStreamingFinished, i, seeder, params, &outputGAM, &outputJSON, &outputGAF, &outputCorrected, &outputCorrectedClipped, &deallocAlns, &stats]() { runComponentMappings(alignmentGraph, readFastqsQueue, readStreamingFinished, i, seeder, params, outputGAM, outputJSON, outputGAF, outputCorrected, outputCorrectedClipped, deallocAlns, stats); });
	}

	for (size_t i = 0; i < params.numThreads; i++)
	{
		threads[i].join();
	}
	assertSetNoRead("Postprocessing");

	allThreadsDone = true;

	GAMwriterThread.join();
	GAFwriterThread.join();
	JSONwriterThread.join();
	correctedWriterThread.join();
	correctedClippedWriterThread.join();
	fastqThread.join();

	if (mummerseeder != nullptr) delete mummerseeder;
	if (minimizerseeder != nullptr) delete minimizerseeder;

	std::string* dealloc;
	while (deallocAlns.try_dequeue(dealloc))
	{
		delete dealloc;
	}

	std::cout << "Alignment finished" << std::endl;
	std::cout << "Input reads: " << stats.reads << " (" << stats.bpInReads << "bp)" << std::endl;
	std::cout << "Seeds found: " << stats.seedsFound << std::endl;
	std::cout << "Seeds extended: " << stats.seedsExtended << std::endl;
	std::cout << "Reads with a seed: " << stats.readsWithASeed << " (" << stats.bpInReadsWithASeed << "bp)" << std::endl;
	std::cout << "Reads with an alignment: " << stats.readsWithAnAlignment << std::endl;
	std::cout << "Alignments: " << stats.alignments << " (" << stats.bpInAlignments << "bp)";
	if (stats.allAlignmentsCount > stats.alignments) std::cout << " (" << (stats.allAlignmentsCount - stats.alignments) << " additional alignments discarded)";
	std::cout << std::endl;
	std::cout << "End-to-end alignments: " << stats.fullLengthAlignments << " (" << stats.bpInFullAlignments << "bp)" << std::endl;
	if (stats.assertionBroke)
	{
		std::cout << "Alignment broke with some reads. Look at stderr output." << std::endl;
	}
}
