# GraphChainer

Here we will version the main project consisting on building a sequence to graph aligner based on Co-linear Chaining.

#### Compilation

To compile, run these:

- Install miniconda https://conda.io/projects/conda/en/latest/user-guide/install/index.html
- `git submodule update --init --recursive`
- `conda env create -f CondaEnvironment.yml`
- `source activate GraphChainer`
- `make bin/GraphChainer`

### Running

Quickstart: `./bin/GraphChainer -t 4 -f reads.fastq -g graph.gfa -a out.gam`

Parameters inherited from GraphAligner
- `-t` Number of threads.
- `-f` Input reads. Format .fasta / .fastq / .fasta.gz / .fastq.gz. You can input multiple files with `-f file1 -f file2 ...` or `-f file1 file2 ...`.
- `-g` Input graph. Format .gfa / .vg.
- `-a` Output file name. Format .gam or .json.

Parameters related to colinear chaining
- `--speed <int>` Default 1. Use 2 or 3 (or larger values) if you want GraphChainer to be faster, but sligthly less accurate.
- `--colinear-split-len <int>` Default 35. The length of the fragments in which the long read is split to create anchors.
- `--colinear-split-gap <int>` Default 35. The distance between consecutive fragments. If `--speed` is set, then always `--colinear-split-gap = --speed * --colinear-split-len`.
- `--colinear-gap 1000` Default 10000. When converting an optimal chain of anchors into an alignment path, split the path if the distance between consecutive anchors is greater than this value.