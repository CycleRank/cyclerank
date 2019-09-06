CycleRank
---------

This repository is a companion for the paper "CycleRank, or There and Back Again: Personalized Relevance Scores from Cyclic Paths on Directed Graphs". 


## CycleRank

`cyclerank.cpp` is the source for the CycleRank algorithm.

It can be compiled with the compilation command:
```
/usr/bin/g++ -Wall -std=c++11 -O2 -flto -pipe -static -s -pthread -o cyclerank cyclerank.cpp
```

### Usage
```
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

## Personalized PageRank


`cyclerank.cpp` is the source for the CycleRank algorithm.

It can be compiled with the compilation command:
```
/usr/bin/g++ pagerank.cpp -Wall -std=c++11 -O2 -static -s -pthread -I/usr/local/include/igraph -L/usr/local/lib -ligraph -o pagerank                                                   1 ↵
```

You want to launch it with the `-w, --wholenetwork` flag to compute it on the whole graph. Furthermore, if you want to compute CheiRank use the `-t, --transposed` flag.

### Usage
```
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