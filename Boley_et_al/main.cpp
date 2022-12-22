/***** A polynomial-delay algorithm for the connector enumeration problem *****/

#define DEBUG
#undef DEBUG

#define MEM
#undef MEM


/***** includes *****/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// C++ libraries
#include <bitset>
#include <chrono>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <vector>
#include <map>
using namespace std;

// custom libraries
#ifdef MEM
#include "mem.cpp"
#endif

#include "define.h"
#include "common.h"
#include "mylib.h"
#include "../statistic/p_value.hpp"


double cpu_time();

/***** declarations *****/
void listGraphKey1(Param P, Tool T, Graph G, BFSTool B, itemset &Item, vector<Solution> &store, Stat stat);
void listGraphBFS(Param P, Tool T, Graph G, BFSTool B,Solution  S, BanList Ban, bitset<ITEM_SIZE> &Item, int key, vector<Solution> &store, Stat stat);
void listGraphDFS(Param P, Tool T, Graph G, BFSTool B,Solution  S, BanList Ban, itemset &Item, int key, vector<Solution> &store, Stat stat);
void nextProcess(Param P, Tool T, Graph G, BFSTool B,int i, bitset<ITEM_SIZE> Item, int key, vector<OwnStack> Q);
void initBFSTool(BFSTool B, int n);


#ifdef DEBUG
int *Dist;


int *F;
#endif

int numAnswer = 0;
int adjTimes = 0;
int k_p = 1;
double checkTime = 0.0;
bool timeout = false; 

/***** main function *****/
int main(int argc, char *argv[]){
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
  readPattern(P, T, G);
  readGraph(P, T, G);
  readPhenotype(P, T);
  readPopulation(P, T);
computeStatisticVal(T);
  initTool(T, G);
  initBFSTool(B, G->n);

  vector<Solution> store;
  _Stat st = _Stat(T, G);
  Stat stat = &st;



  printf("vertices:\t%d\n", G->n);
  printf("edges:\t%d\n", G->m);
  printf("items:\t%d\n", G->items);
  printf("item_density:\t%g\n", G->item_density);

//  printf("%d test Adj\n", G->V[0]->A->v->id);
#ifdef MEM
  printf("BeforeVmSize: %d\nBeforeVmRSS: %d\n", 
	 getVirtualMem(), getPhysicalMem());
#endif
  
#ifdef DEBUG
  //outputGraph(T, G);
  Dist = new int[G->n+1];
  for(int i=0;i<=G->n;i++)
    Dist[i] = 0;


  F = new int[G->n+1];
  for(int i=0;i<=G->n;i++)
    F[i] = 0;
#endif
  itemset Item;
  Item.set();
  printf("run \n");
  /*** run algorithm ***/
  start_time = cpu_time();
  checkTime = start_time;
  auto start = std::chrono::system_clock::now();
  listGraphKey1(P, T, G, B, Item, store, stat);
  cout << "Graphs and p-value" << endl;
    for(auto itr = store.begin(); itr != store.end(); ++itr){
  	  if(stat->p_value(&(*itr))>stat->inverse_threshold(P->alpha, k_p)){
  		  printGraph(&(*itr), T);
  		  cout << "p-value" << Pcmh(P, T, G, &(*itr)) <<endl;
  	  }
    }
  fin_time = cpu_time();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  
#ifdef MEM
  printf("VIRTUAL: %d\nPHYSICAL: %d\n", 
	 getVirtualMem(), getPhysicalMem());
#endif 

#ifdef DEBUG
  for(int i=0;i<=G->n;i++)
    if(Dist[i])
      printf("DIST\t%d\t%d\n",i,Dist[i]);
  delete[] Dist;


  for(int i=0;i<=G->n;i++)
    if(F[i]>1)
      printf("<%d>\n",i);
#endif

  if(timeout){
     cout << "ABORT_BY_TIMELIMIT"<<endl;
  }
  else{
    cout << "TERMINATE_SAFELY"<< endl;
  }

  cout << "all_subgraphs:\t" << numAnswer << "\n";
  cout << "cpu_time:\t" << fin_time-start_time << "\n";
  cout << "elapsed_time:\t" << elapsed_seconds.count() << "\n";
  if(P->reduce){
    printf("reduced_vertices:\t%d\n", T->reduced_vertices);
    printf("reduced_edges:\t%d\n", T->reduced_edges);
  }
  
  /*** postprocess ***/
  delete P;
  delete T;
  delete G;
  return EXIT_SUCCESS;
}

// main algorithm when key=1(i.e. excuting first)
void listGraphKey1(Param P, Tool T, Graph G, BFSTool B, itemset &Item, vector<Solution> &store, Stat stat){
	vector<pair<Solution,BanList> > Q;
	Solution S, bufS;
	BanList Ban;
	itemset nextItem;
	for (int i = 0; i < G->n; ++i) {
		S = bufS;
		if (G->V[i] == NULL)
			continue;
		nextItem.reset();
		int checkNext = getClosure(P, T, G, B, &S, &Ban, G->V[i], Item, &nextItem); // cl(S) でSを更新
		if (checkNext != 0) {
			if (checkNext == 1) {
				Q.push_back(make_pair(S, Ban));
			}
			if (stat->minimal_p_value(&S) > stat->inverse_threshold(P->alpha, k_p)){
				store.push_back(S);

				if (store.size()> k_p){
					k_p++;
					vector<Solution> storebuf;
					auto end_itr = store.end();
					numAnswer = 0;
					for(auto s_itr = store.begin(); s_itr != end_itr; ++s_itr){
						if (stat->minimal_p_value(&(*s_itr)) > stat->inverse_threshold(P->alpha, k_p)){
							storebuf.push_back(*s_itr);
							numAnswer++;
						}
					}
					store = storebuf;
				}
			}
		}
		Ban.push_back(i);
	}

	auto end_itr = Q.end();
	for(auto s_itr = Q.begin(); s_itr != end_itr; ++s_itr){
		if (stat->envelope(&S) > stat->inverse_threshold(P->alpha, k_p)){
			if (P->turnWidth <= 1)
				listGraphDFS(P, T, G, B, s_itr->first,s_itr->second , Item, 2, store, stat);
			else
				listGraphBFS(P, T, G, B, s_itr->first,s_itr->second,Item, 2, store, stat);
		}
	}

}

void listGraphBFS(Param P, Tool T, Graph G, BFSTool B,Solution  S, BanList Ban, itemset &Item, int key, vector<Solution> &store, Stat stat){
	Solution bufS = S;
	vector<pair<Solution,BanList> > Q;
	vector<int> adjVList = adjList(&Ban, G, &S);
	vector<Solution> bufStore;
	auto end = adjVList.end();
	itemset nextItem;
	for (auto itr = adjVList.begin(); itr != end; ++itr){
		S = bufS; //更新されたSを初期に戻す
		int i = *itr;
		int checkNext = getClosure(P, T, G, B, &S, &Ban, G->V[i], Item, &nextItem); // Sを更新
		if (checkNext != 0) {
			if (checkNext == 1)
				Q.push_back(make_pair(S, Ban));
			if (stat->minimal_p_value(&S) < stat->inverse_threshold(P->alpha, k_p)){
				store.push_back(S);
				if (store.size()> k_p){
					k_p++;
					bufStore.clear();
					auto end_itr = store.end();
					numAnswer = 0;
					for(auto s_itr = store.begin(); s_itr != end_itr; ++s_itr){
						if (stat->minimal_p_value(&(*s_itr)) > stat->inverse_threshold(P->alpha, k_p)){
							bufStore.push_back(*s_itr);
							numAnswer++;
						}
					}
					store = bufStore;
				}
			}
		}
		Ban.push_back(i);
	}
	auto end_itr = Q.end();
		for(auto s_itr = Q.begin(); s_itr != end_itr; ++s_itr){
			if (stat->envelope(&S) > stat->inverse_threshold(P->alpha, k_p)){
				if (key == P->turnWidth)
					listGraphDFS(P, T, G, B, (*s_itr).first, (*s_itr).second, Item, key+1, store, stat);
				else
					listGraphBFS(P, T, G, B, (*s_itr).first, (*s_itr).second, Item, key+1, store, stat);
			}
		}
}


void listGraphDFS(Param P, Tool T, Graph G, BFSTool B,Solution  S, BanList Ban, itemset &Item, int key, vector<Solution> &store, Stat stat){
	Solution bufS;
	bufS= S;
	vector<int> adjVList = adjList(&Ban, G, &S);
	vector<Solution> bufStore;
	itemset nextItem;
	auto end = adjVList.end();
	for (auto itr = adjVList.begin(); itr != end; ++itr){
		S = bufS; //更新されたSを初期に戻す
		int i = *itr;
		int checkNext = getClosure(P, T, G, B, &S, &Ban, G->V[i], Item, &nextItem); // Sを更新
		if (checkNext != 0) {
			if (stat->minimal_p_value(&S) > stat->inverse_threshold(P->alpha, k_p)){
				store.push_back(S);
				if (store.size()> k_p){
					k_p++;
					bufStore.clear();
					auto end_itr = store.end();
					numAnswer = 0;
					for(auto s_itr = store.begin(); s_itr != end_itr; ++s_itr){
						if (stat->minimal_p_value(&(*s_itr)) >stat->inverse_threshold(P->alpha, k_p)){
							bufStore.push_back(*s_itr);
							numAnswer++;
						}
					}
					store = bufStore;
				}
			}
			if (checkNext == 1 && stat->envelope(&S) > stat->inverse_threshold(P->alpha, k_p))
				listGraphDFS(P, T, G, B, S, Ban, Item, key+1, store, stat);
		}
		Ban.push_back(i);
	}
}

void initBFSTool(BFSTool B, int n){
  B->marker = 1;
  B->marker_max = MARKER_MAX;
  for(int i=0;i<n;i++){
    B->mark.push_back(0);
    B->bfs_color.push_back(BFS_WHITE);
  }
}
