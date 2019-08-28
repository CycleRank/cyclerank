#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <queue>
#include <stack>
#include <list>
#include <map>
#include <climits>
#include <stdlib.h>     /* exit, EXIT_FAILURE */

extern "C" {
   #include <igraph.h>
}

#include "cxxopts.hpp"
#include "spdlog/spdlog.h"

using namespace std;
namespace spd = spdlog;
namespace opts = cxxopts;

// ***************************************************************************
struct nodo{
  bool active;
  bool blocked;
  int dist;
  vector<int> adj;
  list<int> B;

  nodo(){
    active = true;
    blocked = false;
    dist = -1;
    B.clear();
  }
};


// *************************************************************************
// global variables
shared_ptr<spd::logger> console;

// *************************************************************************
// helper functions
void print_circuit(stack<int> s, vector<int>& new2old) {
  vector<int> tmp;

  printf("--> cycle: ");
  while (!s.empty()) {
    int el = s.top();
    int oldel = new2old[el];
    tmp.push_back(oldel);
    s.pop();
  }

  for (int k=tmp.size()-1; k>=0; k--) {
    printf("%d-", tmp[k]);
  }
  int last = tmp[tmp.size()-1];
  printf("%d\n", last);

}

void write_circuit(stack<int> s, vector<int>& new2old, ofstream& out) {
  vector<int> tmp;

  while (!s.empty()) {
    int el = s.top();
    int oldel = new2old[el];
    tmp.push_back(oldel);
    s.pop();
  }

  for (unsigned k=tmp.size()-1; k>0; k--) {
    out << tmp[k] << " ";
  }
  out << tmp[0] << endl;
}

// count the number of parameters in the first line of the file
// https://stackoverflow.com/a/34665370/2377454
int count_parameters(ifstream& in) {

  int tmp_n;
  int firstline_count =0;
  string line;

  // read first line of file
  getline(in, line);

  stringstream ss(line);

  // read numbers from the line
  while (ss >> tmp_n) {
    ++firstline_count;
  }

  console->debug("firstline_count: {}", firstline_count);

  return firstline_count;
}

// get remapped node
int get_remapped_node_or_fail(int s, map<int,int>& map_old2new ) {
  int newS = -1;

  if ( map_old2new.find(s) == map_old2new.end() ) {
    // Key s not found
    cerr << "Key " << s << " not found in map" << endl;
    exit(EXIT_FAILURE);
  } else {
    // Key s found
    newS = map_old2new[s];
    return newS;
  }

}
// ********** end: helper functions


// *************************************************************************
void transpose_graph(vector<nodo>& g, vector<nodo>& gT) {
  gT.resize(g.size());

  for(unsigned int i=0; i<g.size(); i++) {
    for (int v: g[i].adj) {
      gT[v].adj.push_back(i);
    }
  }
}


void destroy_nodes(vector<nodo>& g, vector<bool>& destroy) {

  for(unsigned int i=0; i<g.size(); i++) {
    if (destroy[i]) {
      g[i].active = false;
      g[i].adj.clear();
    } else {
      for (unsigned int j=0; j<g[i].adj.size(); j++) {
        int v = g[i].adj[j];
        if (destroy[v]) {
          g[i].adj.erase(g[i].adj.begin()+j);
        }
      }
    }
  }
}


void bfs(int source, unsigned int K, vector<nodo>& g) {

  if (!g[source].active) {
    return;
  }

  for (nodo& n: g) {
    n.dist = -1;
  }

  g[source].dist = 0;

  queue<int> q;
  q.push(source);
  int cur;
  while (!q.empty()){
    cur = q.front();
    q.pop();

    // if g[cur].dist == K-1 we can stop
    if(g[cur].dist > (int) K-2) {
      continue;
    }

    for (int v: g[cur].adj) {
      if ((g[v].dist==-1) and (g[v].active)) {

        // neighbor not yet visited, set distance
        g[v].dist = g[cur].dist + 1;
        q.push(v);
      }
    }
  }
}


int main(int argc, const char* argv[]) {

  // *************************************************************************
  // initialize logger
  try {
    console = spd::stdout_color_mt("console");
  }
  // exceptions thrown upon failed logger init
  catch (const spd::spdlog_ex& ex) {
    cerr << "Log init failed: " << ex.what() << endl;
    return 1;
  }
  // ********** end: logger

  // *************************************************************************
  // parse command-line options
  opts::Options* options;
  string input_file="input.txt";
  string output_file="output.txt";
  double alpha=0.85;
  int cliS = -1;
  int cliK = -1;
  bool verbose = false;
  bool debug = false;
  bool help = false;
  bool transposed = false;
  bool undirected = false;
  bool directed = true;
  bool wholenetwork = false;
  bool forcebfstransposed = false;

  try {
    options = new cxxopts::Options(argv[0]);

    options->add_options()
      ("a,alpha", "Damping factor (alpha).",
       cxxopts::value<double>(alpha),
       "ALPHA"
       )
      ("f,file", "Input file.",
       cxxopts::value<string>(input_file),
       "FILE"
       )
      ("v,verbose", "Enable logging at verbose level.",
       cxxopts::value(verbose))
      ("d,debug", "Enable logging at debug level.",
       cxxopts::value(debug))
      ("h,help", "Show help message and exit.",
       cxxopts::value(help))
      ("k,maxloop", "Set max loop length (K).",
       cxxopts::value(cliK),
       "K"
       )
      ("o,output", "Output file.",
       cxxopts::value<string>(output_file),
       "OUTPUT_FILE"
       )
      ("s,source", "Set source node (S).",
       cxxopts::value(cliS),
       "S"
       )
      ("t,transposed", "Run on the transposed network (incompatible with -u).",
       cxxopts::value(transposed))
      ("u,undirected", "Run on the undirected network (incompatible with -t).",
       cxxopts::value(undirected))
      ("b,force-bfs-transposed", "Force running the second BFS on the " \
       "transposed network. This is needed only if -u is specified, " \
       "otherwise it is effectively ignored.",
       cxxopts::value(forcebfstransposed))
      ("w,whole-network", "Run on the whole network (ignore K).",
       cxxopts::value(wholenetwork))
      ;

    auto arguments = options->parse(argc, argv);
  } catch (const cxxopts::OptionException& e) {
    cerr << "Error parsing options: " << e.what() << endl;
    exit (EXIT_FAILURE);
  }

  // if help option is activated, print help and exit.
  if(help) {
    cout << options->help({""}) << endl;
    exit(0);
  }

  if(transposed && undirected) {
    cerr << "Error: options -t (transposed) and -u (undirected) are " \
         << "mutually exclusive." << endl;
    exit (EXIT_FAILURE);
  }

  if(alpha <= 0) {
    cerr << "Error: the damping factor specified with -a (alpha) must be" \
         << "positive." << endl;
    exit (EXIT_FAILURE);
  }

  // use a variable called directed to indicate if the graph is directed
  directed = !undirected;
  // ********** end: parse command-line options

  // *************************************************************************
  // start logging
  // set logging level based on option from CLI
  if (debug) {
    spd::set_level(spd::level::debug);
  } else if (verbose) {
    spd::set_level(spd::level::info);
  } else {
    spd::set_level(spd::level::warn);
  }

  console->info("Log start!");
  console->debug("input_file: {}", input_file);
  console->debug("verbose: {}", verbose);
  console->debug("debug: {}", debug);
  console->debug("alpha: {}", alpha);
  console->debug("transposed: {}", transposed);
  console->debug("undirected: {}", undirected);
  console->debug("directed: {}", directed);
  console->debug("whole-network: {}", wholenetwork);
  // ********** end: start logging

  // *************************************************************************
  // start algorithm
  int S = -1, newS = -1;
  unsigned int N = 0, M = 0, K = 0;
  vector<nodo> grafo;

  map<int,int> old2new;
  vector<int> new2old;

  int count_destroied = 0;

  // *************************************************************************
  // read input
  {
    ifstream in(input_file);
    int tmpS = -1;
    int tmpK = -1;
    int nparam = 0;

    if(in.fail()){
      cerr << "Error! Could not open file: " << input_file << endl;
      exit(EXIT_FAILURE);
    }

    nparam = count_parameters(in);
    in.close();

    in.open(input_file);
    if(nparam == 4) {
      in >> N >> M >> tmpS >> tmpK;
    } else if(nparam == 2) {
      in >> N >> M;
    } else {
      cerr << "Error! Error while reading file (" << input_file \
          << "), unexpected number of parameters" << endl;
      exit(EXIT_FAILURE);
    }

    if(cliS == -1) {
      S = tmpS;
    } else {
      S = cliS;
    }

    if(cliK == -1) {
      K = (unsigned int) tmpK;
    } else {
      K = (unsigned int) cliK;
    }

    assert( (N > 0 && M > 0) \
            && "N and M must be positive." );

    assert( (K > 0 && S >= 0) \
            && "K must be positive and S must be non-negative." );

    console->info("N: {}", N);
    console->info("M: {}", M);
    console->info("S: {}", S);
    console->info("K: {}", K);

    console->debug("reading graph...");
    grafo.resize(N);
    for(unsigned int j=0; j<M; j++) {
      int s, t;
      in >> s >> t;

      // check that we are not inserting duplicates
      if (find(grafo[s].adj.begin(), \
               grafo[s].adj.end(), \
               t) == grafo[s].adj.end()) {
        grafo[s].adj.push_back(t);

        if(!directed) {
          // we still need to check for duplicates since we are starting
          // from a directed network and if we a double link, i.e.:
          //   1 -> 2
          //   2 -> 1
          // this would become a duplicate when considering the undirected
          // version
          if (find(grafo[t].adj.begin(), \
                   grafo[t].adj.end(), \
                   s) == grafo[t].adj.end()) {
            grafo[t].adj.push_back(s);
          }
        }
      }
    }
    console->debug("--> read graph");
    in.close();
  }
  // ********** end: read input

  if(!wholenetwork) {
    // ***********************************************************************
    // Step 1: BFS on g
    {
      console->info("Step 1. BFS");
      vector<bool> destroy(N, false);

      bfs(S, K, grafo);

      for(unsigned int i=0; i<N; i++) {
        // se il nodo si trova distanza maggiore di K-1 o non raggiungibile
        if((grafo[i].dist == -1) or (grafo[i].dist > (int) K-1)) {
          destroy[i] = true;
          count_destroied++;
        }
      }

      int remaining = N-count_destroied;
      console->info("nodes: {}", N);
      console->info("destroyed: {}", count_destroied);
      console->info("remaining: {}", remaining);
      new2old.resize(remaining);

      int newindex = -1;
      for(unsigned int i=0; i<N; i++) {
        if(!destroy[i]) {
          newindex++;
          new2old[newindex] = i;
          old2new.insert(pair<int,int>(i,newindex));
          /*
          console->debug("old2new.insert(pair<int,int>({0}, {1}))",
                         i,
                         newindex);
          */
        }
      }

      destroy_nodes(grafo, destroy);

      vector<nodo> tmpgrafo;
      tmpgrafo.resize(remaining);

      int newi = -1;
      int newv = -1;

      for(unsigned int i=0; i<N; i++) {
        if(grafo[i].active) {
          newi = old2new[i];
          tmpgrafo[newi].dist = grafo[i].dist;
          tmpgrafo[newi].active = true;
          for (int v: grafo[i].adj) {
            if(grafo[v].active) {
              newv = old2new[v];
              tmpgrafo[newi].adj.push_back(newv);
            }
          }
        }
      }

      grafo.clear();
      grafo.swap(tmpgrafo);

      destroy.clear();
    }
    // ********** end: Step 1

    vector<bool> destroy(grafo.size(), false);


    // ***********************************************************************
    // get remapped source node (S)
    newS = get_remapped_node_or_fail(S, old2new);
    console->info("S: {0}, newS: {1}", S, newS);
    // ********** end: get remapped source node (S)

    // ***********************************************************************
    // Step 2: BFS on g^T
    if(directed || forcebfstransposed) {
      {
        console->info("Step 2.: BFS on g^T");

        if(directed) {
          // graph is directed, this is the regular scenario
          vector<nodo> grafoT;
          transpose_graph(grafo, grafoT);

          bfs(newS, K, grafoT);

          for(unsigned int i=0; i<grafo.size(); i++) {
            if((grafo[i].dist == -1) or (grafoT[i].dist == -1) or \
                (grafo[i].dist + grafoT[i].dist > (int) K)) {
              // console->debug("destroied node: {0:d}\n", i);
              destroy[i] = true;
              count_destroied++;
            }
          }

        } else {
          // graph is undirected, we are here because we are forcing this
          // step

          // we can be smart and instead of running the same visit twice
          // (we are on an undirected network, we transpose it and we run a
          // BFS from the same node), we know that:
          //   \foreach i grafo[i].dist == grafoT[i].dist
          //
          // hence the check over the distance we do above:
          //   (grafo[i].dist == -1) || (grafoT[i].dist == -1) ||
          //       (grafo[i].dist + grafoT[i].dist > K)
          //
          // is effectively simplified to:
          //   (grafo[i].dist == -1) || (2*grafo[i].dist > K)
          //
          for(unsigned int i=0; i<grafo.size(); i++) {
            if((grafo[i].dist == -1) or (2*grafo[i].dist > (int) K)) {
              // console->debug("destroied node: {0:d}\n", i);
              destroy[i] = true;
              count_destroied++;
            }
          }

        }

        int remaining = N-count_destroied;

        console->info("nodes: {}", N);
        console->info("destroyed: {}", count_destroied);
        console->info("remaining: {}", remaining);

        destroy_nodes(grafo, destroy);
      }
      // ********** end: Step 2

      int remaining = N-count_destroied;
      map<int,int> tmp_old2new;
      vector<int> tmp_new2old;
      tmp_new2old.resize(remaining);

      {
        int newindex = -1;
        int oldi = -1;
        for(unsigned int i=0; i<grafo.size(); i++) {
          if(!destroy[i]) {
            newindex++;

            oldi = new2old[i];

            tmp_new2old[newindex] = oldi;
            tmp_old2new.insert(pair<int,int>(oldi, newindex));
          }
        }
      }

      {
        int oldi = -1, tmpnewi = -1;
        int oldv = -1, tmpnewv = -1;

        vector<nodo> tmpgrafo;
        tmpgrafo.resize(remaining);
        for(unsigned int i=0; i<grafo.size(); i++) {
          if(grafo[i].active) {

            oldi = new2old[i];
            tmpnewi = tmp_old2new[oldi];
            tmpgrafo[tmpnewi].dist = grafo[i].dist;
            tmpgrafo[tmpnewi].active = true;
            for (int v: grafo[i].adj) {
              if(grafo[v].active) {
                oldv = new2old[v];
                tmpnewv = tmp_old2new[oldv];
                tmpgrafo[tmpnewi].adj.push_back(tmpnewv);
              }
            }
          }
        }

        grafo.clear();
        grafo.swap(tmpgrafo);
        destroy.clear();
      }

      new2old.clear();
      old2new.clear();

      tmp_new2old.swap(new2old);
      tmp_old2new.swap(old2new);

      // *********************************************************************
      // get remapped source node (S)
      newS = get_remapped_node_or_fail(S, old2new);
      console->info("S: {0}, newS: {1}", S, newS);
      // ********** end: get remapped source node (S)
    }
  } else {
    console->info("Running on the whole network");
    // we have the original network
    newS = S;
  }
  /* *************************************************************************
  * Personalized pagerank from igraph
  *   http://igraph.org/c/doc/
  *       igraph-Structural.html#igraph_personalized_pagerank
  *
  * int igraph_personalized_pagerank(
  *    const igraph_t *graph, 
  *    igraph_pagerank_algo_t algo, igraph_vector_t *vector,
  *    igraph_real_t *value, const igraph_vs_t vids,
  *    igraph_bool_t directed, igraph_real_t damping, 
  *    igraph_vector_t *reset,
  *    const igraph_vector_t *weights,
  *    void *options);
  * *************************************************************************/
  console->info("Pagerank (alpha={})", alpha);

  igraph_t igrafo;
  igraph_vector_t iedges;
  vector<nodo> grafoT;
  vector<nodo> grafoU;
  igraph_real_t pr_alpha(alpha);

  console->debug("Calculating the Pagerank on the graph: ");
  if(!transposed) {
    console->debug("  * on the input graph");
  } else {
    console->debug("  * on the transposed graph");
    transpose_graph(grafo, grafoT);
    grafo.swap(grafoT);

    // deallocate vector (of nodes)
    grafoT.clear();
  }

  if(directed) {
    console->debug("  * on the directed graph");
  } else {
    console->debug("  * on the undirected graph");
  }

  // count number of edges
  unsigned int num_edges = 0;
  for(unsigned int i=0; i<grafo.size(); i++) {
    num_edges += grafo[i].adj.size();
  }

  console->debug("num_edges: {0}", num_edges);

  // initialize vertices vector 2*num_edges
  igraph_vector_init(&iedges, 2*num_edges);

  int ec = 0;
  for(unsigned int i=0; i<grafo.size(); i++) {
    for (int v: grafo[i].adj) {
      VECTOR(iedges)[ec]=i;
      VECTOR(iedges)[ec+1]=v;

      ec = ec + 2;
    }
  }

  // we get the size of the graph and the we clear
  unsigned int num_nodes = grafo.size();
  grafo.clear();

  // int igraph_create(igraph_t *graph, const igraph_vector_t *edges,
  //   igraph_integer_t n, igraph_bool_t directed);
  igraph_create(&igrafo, &iedges, num_nodes, directed);

  igraph_vector_t pprscore, reset;

  // init result vector
  igraph_vector_init(&pprscore, 0);

  // reset vector
  igraph_vector_init(&reset, num_nodes);
  igraph_vector_fill(&reset, 0);

  /*
  * int igraph_personalized_pagerank(
  *    const igraph_t *graph,
  *    igraph_pagerank_algo_t algo, igraph_vector_t *vector,
  *    igraph_real_t *value, const igraph_vs_t vids,
  *    igraph_bool_t directed, igraph_real_t damping,
  *    igraph_vector_t *reset,
  *    const igraph_vector_t *weights,
  *    void *options);
  *
  * algo: IGRAPH_PAGERANK_ALGO_PRPACK is the recommended implementation
  *       http://igraph.org/c/doc/igraph-Structural.html#igraph_pagerank_algo_t
  */

  // jump probability to (remapped) S is 1
  VECTOR(reset)[newS]=1.0;

  int ret = -1;
  ret=igraph_personalized_pagerank(
     &igrafo,                         // const igraph_t *graph
     IGRAPH_PAGERANK_ALGO_PRPACK,     // igraph_pagerank_algo_t algo
     &pprscore,                       // igraph_vector_t *vector
     0,                               // igraph_real_t *value
     igraph_vss_all(),                // const igraph_vs_t vids
     directed,                        // igraph_bool_t directed
     pr_alpha,                        // igraph_real_t damping
     &reset,                          // igraph_vector_t *reset
     0,                               // const igraph_vector_t *weights,
     0                                // void *options
     );
  console->debug("SSPPR ret: {0:d}", ret);

  igraph_vector_destroy(&reset);
  igraph_destroy(&igrafo);

  FILE* outfp;
  outfp = fopen(output_file.c_str(), "w+");
  for (unsigned int i=0; i<num_nodes; i++) {
    int oldi;
    if(!wholenetwork) {
      oldi = new2old[i];
    } else {
      oldi = i;
    }
    fprintf(outfp, "score(%d):\t%.10f\n", oldi, VECTOR(pprscore)[i]);
  }

  console->info("Log stop!");
  exit (EXIT_SUCCESS);
}
