# CycleRank

This repository is a companion for the paper "CycleRank, or There and Back
Again: Personalized Relevance Scores from Cyclic Paths on Directed Graphs".

## Dependencies

The code needs the following libraries as dependencies:

* [igraph][igraph] v0.7.1
* [spdlog][spdlog] v0.17.0
* [cxxopts][cxxopts] v2.1.0

## Assumptions

### Input Format

The input is supposed to be M+1 lines long.

```plaintext
N M [S K]

```

* the first line can be 2 or 4 integer numbers:

  ```plaintext
  N M [S K]
  s_1 t_1
  ...
  s_M t_M
  ```

  where:
  * `N` is the total number of nodes
  * `M` is the total number of edges
  * `S` is the id of the seed node
  * `K` is the maximum length of cycles to consider

* the following M lines are two integer:

```plaintext
s_i t_i
```

  where `s_i`, `t_i` are the ids of the source and target nodes
  and `0 <= s_i,t_i <= N`

### Input Characteristics

The code assumes that node indexes are consecutive, starting from 0.
In other words, the input must be relabeled so that the there are no gaps in indexes.

## Algorithms

### CycleRank

`cyclerank.cpp` is the source for the CycleRank algorithm.

It can be compiled with the compilation command:

```bash
/usr/bin/g++ -Wall -std=c++11 -O2 -flto -pipe -static -s -pthread \
  -o cyclerank cyclerank.cpp
```

#### Usage

```bash
Usage:
  ./cyclerank [OPTION...]

  -f, --file FILE           Input file.
  -v, --verbose             Enable logging at verbose level.
  -d, --debug               Enable logging at debug level.
  -h, --help                Show help message and exit.
  -k, --maxloop K           Set max loop length (K).
  -o, --output OUTPUT_FILE  Output file.
  -s, --source S            Set source node (S).
```

### Personalized PageRank

`cyclerank.cpp` is the source for the CycleRank algorithm.

It can be compiled with the compilation command:

```bash
/usr/bin/g++ pagerank.cpp -Wall -std=c++11 -O2 -static -s -pthread \
  -I/usr/local/include/igraph -L/usr/local/lib -ligraph -o pagerank
```

You want to launch it with the `-w, --wholenetwork` flag to compute it on the
whole graph. Furthermore, if you want to compute CheiRank use the `-t,
--transposed` flag.

#### Usage

```bash
Usage:
  ./pagerank [OPTION...]

  -a, --alpha ALPHA           Damping factor (alpha).
  -f, --file FILE             Input file.
  -v, --verbose               Enable logging at verbose level.
  -d, --debug                 Enable logging at debug level.
  -h, --help                  Show help message and exit.
  -k, --maxloop K             Set max loop length (K).
  -o, --output OUTPUT_FILE    Output file.
  -s, --source S              Set source node (S).
  -t, --transposed            Run on the transposed network (incompatible
                              with -u).
  -u, --undirected            Run on the undirected network (incompatible
                              with -t).
  -b, --force-bfs-transposed  Force running the second BFS on the transposed
                              network. This is needed only if -u is specified,
                              otherwise it is effectively ignored.
  -w, --whole-network         Run on the whole network (ignore K).
```

[igraph]: https://github.com/igraph/igraph
[spdlog]: https://github.com/gabime/spdlog
[cxxopts]: https://github.com/jarro2783/cxxopts
