/***** A polynomial-delay algorithm for the connector enumeration problem *****/

#define DEBUG
#undef DEBUG

#define MEM
#undef MEM

/***** includes *****/
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// C++ libraries
#include <bitset>
#include <chrono>
#include <iostream>
#include <map>
#include <queue>
#include <unordered_map>
#include <vector>
#include <tuple>
using namespace std;

// custom libraries
#ifdef MEM
#include "mem.cpp"
#endif

#include "../statistic/p_value.hpp"
#include "common.h"
#include "define.h"
#include "mylib.h"

double cpu_time();

/***** declarations *****/
void listGraphKey1(Param P, Tool T, Graph G, BFSTool B, itemset &Item,
                   multimap<double, Solution> &container, Stat stat, int depth);
void listGraphBFS(Param P, Tool T, Graph G, BFSTool B, Solution S, BanList Ban,
                  bitset<ITEM_SIZE> &Item, int key, multimap<double, Solution> &container,
                  Stat stat, int depth);
void listGraphDFS(Param P, Tool T, Graph G, BFSTool B, Solution S, BanList Ban,
                  itemset &Item, int key, multimap<double, Solution> &container, Stat stat, int depth);
void updateStore(Param P, Solution S, multimap<double, Solution> &container, Stat stat);
void findSignificants(Param P, multimap<double, Solution> &container, Stat stat);
void nextRecursive(Param P, Tool T, Graph G, BFSTool B, itemset &Item, int key,
                   multimap<double, Solution> &container, Stat stat,
                   vector<tuple<Solution, BanList, itemset>> &Q, int depth);
void initBFSTool(BFSTool B, int n);

#ifdef DEBUG
int *Dist;

int *F;
#endif

int numAnswer = 0;
int numSignificant = 0;
int numNode = 0;
int adjTimes = 0;
int k_p = 1;
double checkTime = 0.0;

bool timeout = false;
std::chrono::system_clock::time_point start_t;
std::chrono::system_clock::time_point end_t;


/***** main function *****/
int main(int argc, char *argv[]) {
  Param P;
  Tool T;
  Graph G;
  BFSTool B;
  double start_time, fin_time;

  /*** read the arguments ***/
  checkArgs(argc, argv);
  P = new struct _Param;
  T = new struct _Tool;
  G = new struct _Graph;
  B = new struct _BFSTool;

  readArgs(P, G, argc, argv);
  cpu_time();

  /*** read the pattern file and the graph file ***/
  preprocess(P, T, G);
  //cout << "!" <<endl;
  readPattern(P, T, G);
  //cout << "!!" <<endl;
  readGraph(P, T, G);
  //cout << "!!!" <<endl;
  //readPhenotype(P, T);
  //readPopulation(P, T);
  computeStatisticVal(T);
  initTool(T, G);
  initBFSTool(B, G->n);

  multimap<double, Solution> container;
  _Stat st = _Stat(T, G);
  Stat stat = &st;

  printf("vertices:\t%d\n", G->n);
  printf("edges:\t%d\n", G->m);
  printf("items:\t%d\n", G->items);
  printf("item_density:\t%g\n", G->item_density);

#ifdef MEM
  printf("BeforeVmSize: %d\nBeforeVmRSS: %d\n", getVirtualMem(),
         getPhysicalMem());
#endif

#ifdef DEBUG
  // outputGraph(T, G);
  Dist = new int[G->n + 1];
  for (int i = 0; i <= G->n; i++)
    Dist[i] = 0;

  F = new int[G->n + 1];
  for (int i = 0; i <= G->n; i++)
    F[i] = 0;
#endif
  itemset Item;
  Item.set();
  printf("run \n");
  /*** run algorithm ***/
  start_time = cpu_time();
  checkTime = start_time;
  numAnswer = 0; //container.size();
  start_t = std::chrono::system_clock::now();
  listGraphKey1(P, T, G, B, Item, container, stat, 0);
  //findSignificants(P, container, stat);
  fin_time = cpu_time();
  end_t = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds;
  elapsed_seconds = end_t - start_t;

  cout << "k_p: " << k_p << endl;
  cout << "Significants: " << container.size() << endl;
  cout << "outname: " << P->outname << endl;
  

#ifdef MEM
  printf("VIRTUAL: %d\nPHYSICAL: %d\n", getVirtualMem(), getPhysicalMem());
#endif

#ifdef DEBUG
  for (int i = 0; i <= G->n; i++)
    if (Dist[i])
      printf("DIST\t%d\t%d\n", i, Dist[i]);
  delete[] Dist;

  for (int i = 0; i <= G->n; i++)
    if (F[i] > 1)
      printf("<%d>\n", i);
#endif

  if (timeout) {
    cout << "ABORT_BY_TIMELIMIT" << endl;
  } else {
    cout << "TERMINATE_SAFELY" << endl;
  }

  cout << "all_subgraphs:\t" << numAnswer << "\n";
  cout << "all_ignificants:\t" << numSignificant << "\n";
  cout << "all_nodes\t" << numNode << endl;
  cout << "cpu_time:\t" << fin_time - start_time << "\n";
  cout << "elapsed_time:\t" << elapsed_seconds.count() << "\n";
  if (P->reduce) {
    printf("reduced_vertices:\t%d\n", T->reduced_vertices);
    printf("reduced_edges:\t%d\n", T->reduced_edges);
  }

  writeSignificantsToFile(std::string(P->outname), container);

  /*** postprocess ***/
  delete P;
  delete T;
  delete G;
  return EXIT_SUCCESS;
}

// main algorithm when key=1(i.e. excuting first)
void listGraphKey1(Param P, Tool T, Graph G, BFSTool B, itemset &Item,
                   multimap<double, Solution> &container, Stat stat, int depth) {
  vector<tuple<Solution, BanList, itemset>> Q;
  Solution S, bufS;
  BanList Ban;
  itemset nextItem;
  
  for (int i = 0; i < G->n; ++i) {
    std::chrono::duration<double> elapsed;
    auto now_t = std::chrono::system_clock::now();
    elapsed = now_t-start_t;
    //printf("%g\t%g\t%g\n", elapsed, elapsed.count(),P->tlim);
    if(elapsed.count() > P->tlim){
      timeout = true;
      return;
    }
    
    S = bufS;
    if (G->V[i] == NULL)
      continue;
    nextItem.reset();
    int checkNext = getClosure(P, T, G, B, &S, &Ban, G->V[i], Item,
                               &nextItem); // cl(S) でSを更新
    if (checkNext != 0) {
      if (checkNext == 1) {
        Q.emplace_back(S, Ban, nextItem);
      } else {
        // This node does not have any candidate of child
        numAnswer++;
      }
      updateStore(P, S, container, stat);
    }
    Ban.push_back(i);
  }
  nextRecursive(P, T, G, B, Item, 2, container, stat, Q, depth + 1);
}

void listGraphBFS(Param P, Tool T, Graph G, BFSTool B, Solution S, BanList Ban,
                  itemset &Item, int key, multimap<double, Solution> &container, Stat stat, int depth) {
  Solution bufS = S;
  vector<tuple<Solution, BanList, itemset>> Q;
  vector<int> adjVList = adjList(&Ban, G, &S);
  vector<Solution> bufStore;
  auto end = adjVList.end();
  itemset nextItem;

  for (auto itr = adjVList.begin(); itr != end; ++itr) {
    S = bufS; // 更新されたSを初期に戻す
    int i = *itr;
    int checkNext =
        getClosure(P, T, G, B, &S, &Ban, G->V[i], Item, &nextItem); // Sを更新
    if (checkNext != 0) {
      if (checkNext == 1) {
        Q.emplace_back(S, Ban, nextItem);
      } else {
        // This node does not have any candidate of child
        numAnswer++;
      }
      updateStore(P, S, container, stat);
    }
    Ban.push_back(i);
  }
  auto end_itr = Q.end();
  nextRecursive(P, T, G, B, Item, key, container, stat, Q, depth+1);
}

void listGraphDFS(Param P, Tool T, Graph G, BFSTool B, Solution S, BanList Ban,
                  itemset &Item, int key, multimap<double, Solution> &container, Stat stat, int depth) {
  Solution bufS;
  bufS = S;
  vector<int> adjVList = adjList(&Ban, G, &S);
  vector<Solution> bufStore;
  itemset nextItem;
  auto end = adjVList.end();

  for (auto itr = adjVList.begin(); itr != end; ++itr) {
    S = bufS; // 更新されたSを初期に戻す
    int i = *itr;
    int checkNext =
        getClosure(P, T, G, B, &S, &Ban, G->V[i], Item, &nextItem); // Sを更新
    if (checkNext != 0) {
      if (key % 2 == 0) {
        updateStore(P, S, container, stat);
        numAnswer++;
      }
      if (checkNext == 1 &&
          stat->envelope(&S) > stat->inverse_threshold(P->alpha, k_p)) {
        listGraphDFS(P, T, G, B, S, Ban, Item, key + 1, container, stat, depth+1);
      }
      if (key % 2 != 0) {
        updateStore(P, S, container, stat);
        numAnswer++;
      }
    }
    Ban.push_back(i);
  }
}

void initBFSTool(BFSTool B, int n) {
  B->marker = 1;
  B->marker_max = MARKER_MAX;
  for (int i = 0; i < n; i++) {
    B->mark.push_back(0);
    B->bfs_color.push_back(BFS_WHITE);
  }
}

void updateStore(Param P, Solution S, multimap<double, Solution> &container, Stat stat) {
  // auto value = stat->minimal_p_value(&S);
  // HERE 2
  // auto threshold = stat->inverse_threshold(P->alpha, k_p);
  // if (value > threshold) {
  
  /*container.emplace(value, S);*/
    
    // if (container.size() > k_p) {
    //   k_p++;
    //   auto itr = container.begin();
    //   auto end = container.lower_bound(threshold);
    //   // Remove candidates whose p-value is greater than threshold.
    //   // Note that we use the inverse function of the survival function
    //   // for the importance comparisons, so we are removing those with smaller values.
    //   container.erase(itr, end);
    // }
  // }
}

void findSignificants(Param P, multimap<double, Solution> &container, Stat stat) {
  vector<double> buf;
  multimap<double, Solution> bufContainer;
  double th = stat->inverse_threshold(P->alpha, k_p);
  auto end_itr = container.end();
  for (auto c_itr = container.begin(); c_itr != end_itr; ++c_itr) {
    double p = stat->p_value(&(c_itr->second));
    if (p > th) {
      auto p_value = stat->survival_function(p);
      bufContainer.emplace(p_value, c_itr->second);
    } else {
      buf.push_back(p);
    }
  }
  numSignificant = bufContainer.size();
  cout << "TH: " << th << endl;
  container = bufContainer;
}

void nextRecursive(Param P, Tool T, Graph G, BFSTool B, itemset &Item, int key,
                   multimap<double, Solution> &container, Stat stat,
                   vector<tuple<Solution, BanList, itemset>> &Q, int depth) {
  auto end_itr = Q.end();
  numNode++;


    

  for (auto s_itr = Q.begin(); s_itr != end_itr; ++s_itr) {

    
  /*** process timeout ***/
  std::chrono::duration<double> elapsed;
  auto now_t = std::chrono::system_clock::now();
  elapsed = now_t-start_t;
  if(elapsed.count() > P->tlim){
    timeout = true;
    return;
  }
  
    // HERE 1
    // if (stat->envelope(&(s_itr->first)) >
    //     stat->inverse_threshold(P->alpha, k_p)) {
      auto [S, Ban, nextItem] = *s_itr;
      if (key % 2 == 0) {
        numAnswer++;
      }
      if (key == P->turnWidth) {
        listGraphDFS(P, T, G, B, S, Ban, nextItem, key + 1,
                     container, stat, depth);
      } else {
        listGraphBFS(P, T, G, B, S, Ban, nextItem, key + 1,
                     container, stat, depth);
      }
      if (key % 2 != 0) {
        numAnswer++;
      }
    // }
  }
}
