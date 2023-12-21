/***** A polynomial-delay algorithm for the connector enumeration problem *****/

#define DEBUG
#undef DEBUG

#define MEM
#undef MEM

int timeout=0;

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
using namespace std;

// custom libraries
#ifdef MEM
#include "mem.cpp"
#endif

#include "define.h"
#include "common.h"
#include "mylib.h"
#include "mt19937ar.h"

double cpu_time();

/***** declarations *****/
HugePositive family_tree(Param P, Tool T, Graph G);
HugePositive descendants(Param P, Tool T, Graph G, BFSTool B, Component Parent, int d);
void initBFSTool(BFSTool B, int n);

#ifdef DEBUG
int *Dist;


int *F;
#endif

double start_time, fin_time;

/***** main function *****/
int main(int argc, char *argv[]){
  Param P;
  Tool T;
  Graph G;
  HugePositive numConnectors;
  
  /*** read the arguments ***/
  checkArgs(argc, argv);
  P = new struct _Param;
  T = new struct _Tool;
  G = new struct _Graph;
  readArgs(P, G, argc, argv);
  cpu_time();

  /*** initialize ***/
  init_genrand(P->seed);

  /*** read the pattern file and the graph file ***/
  preprocess(P, T, G);
  readPattern(P, T, G);
  readGraph(P, T, G);
  initTool(T, G);

  printf("vertices:\t%d\n", G->n);
  printf("edges:\t%d\n", G->m);
  printf("edge_density:\t%g\n", getDensity(G->n, G->m));
  printf("items:\t%d\n", G->items);
  printf("item_density:\t%g\n", G->item_density);
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

  
  /*** run algorithm ***/
  start_time = cpu_time();
  auto start = std::chrono::system_clock::now();
  numConnectors = family_tree(P, T, G);
  fin_time = cpu_time();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;

  if(timeout)
    printf("ABORT_BY_TIMELIMIT\n");
  else
    printf("TERMINATE_SAFELY\n");
  
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

  
  cout << "all_subgraphs:\t" << numConnectors << "\n";
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


HugePositive family_tree(Param P, Tool T, Graph G){
  HugePositive numConnectors = 0,desc;
  BFSTool B;
  
  /*** preparation ***/
  B = new struct _BFSTool;
  initBFSTool(B, G->n);

  // connectors with empty item set
  {
    vector<Component> Cmax;
    int numCmax;
    numCmax = getCmax(T, G, B, T->seq_all, Cmax);
    for(int z=0;z<numCmax;z++){
      if(isCommonItemEmpty(T, G, Cmax[z])){
	numConnectors++;
#ifdef DEBUG
	int num_debug;
	num_debug = Cmax[z]->seq.size();
	Dist[num_debug]++;
#endif
      }
    }
    for(int z=0;z<numCmax;z++)
      delete Cmax[z];

    cout << "empty_itemset_subgraphs:\t" << numConnectors << "\n";
  }
  
  // connectors with nonempty item set
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int k=0;k<G->items;k++){
    vector<Component> Cmax;
    int numCmax;
    Graph S = new struct _Graph;

    if(cpu_time()-start_time > P->tlim && P->tlim>0){
      timeout = 1;
      break;
    }

    
#ifdef _OPENMP
    BFSTool B;
    B = new struct _BFSTool;
    initBFSTool(B, G->n);
    if(k==0)
      printf("num_threads:\t%d\n", omp_get_num_threads());
#endif
    getSubgraph(T, G, S, k);
    numCmax = getCmax(T, S, B, T->seq_item[k], Cmax);
    
    for(int z=0;z<numCmax;z++){

      if(isMinItemK(T, G, Cmax[z], k)){

#ifdef _OPENMP
#pragma omp atomic
#endif
	numConnectors++;
	
#ifdef DEBUG
	int num_debug;
	num_debug = Cmax[z]->seq.size();
	Dist[num_debug]++;
#endif

	if(k<G->items-1){
	  getItemIntersection(T, S, Cmax[z]);
	  desc = descendants(P, T, S, B, Cmax[z], 2);
	  
#ifdef _OPENMP
#pragma omp atomic
#endif
	  numConnectors += desc;
	}
      }
    }
    for(int z=0;z<numCmax;z++){
      if(Cmax[z]->I!=NULL)
	delete[] Cmax[z]->I;
      delete Cmax[z];
    }
    delSubgraph(S);
    delete S;
#ifdef _OPENMP
    delete B;
#endif
  }
  delete B;
  return numConnectors;
}


HugePositive descendants(Param P, Tool T, Graph G, BFSTool B, Component Parent, int d){
  bool first = true;
  HugePositive numConnectors = 0;
  vector<Component> Cmax;
  vector<VertexIDSeq> Pseq;
  int numCmax,i;

  for(int j=Parent->minitem+1;j<G->items;j++){
    if(hasItem(T->B, Parent->I, j))
      continue;
    numCmax = getCmaxByItem(T, G, B, Parent, j, Cmax);
    for(int z=0;z<numCmax;z++){
      if(isMinItemK(T, G, Cmax[z], Parent->minitem)==false)
	continue;
      getItemIntersection(T, G, Cmax[z]);
      i = getDiffMinItem(T, G, Parent->minitem, Parent, Cmax[z]);
    
      if(j != i)
	continue;

      /*****/
      if(first){
	bool minitem=true;
	int len=0;
	for(auto itr=Parent->Itv.begin(); itr!=Parent->Itv.end(); itr++){
	  int p = *itr;
	  BitString b = Parent->I[p];
	  int z,l=0;
	  while(1){
	    l = getFirstBit(T->H, T->T, b, l, BIT_LENGTH-1);
	    if(l==Undef)
	      break;
	    z = p*BIT_LENGTH+l;
	    if(minitem){
	      Pseq.push_back(T->seq_item[Parent->minitem]);
	      minitem = false;
	    }
	    else{
	      VertexIDSeq seq;
	      seq = getSubseqByItem(T, G, z, Pseq[len-1]);
	      Pseq.push_back(seq);
	    }	    
	    len++;
	    l++;
	    if(l==BIT_LENGTH)
	      break;
	  }
	}
	first = false;
      }
      /*****/

      if(!isParent(T, G, B, Cmax[z], Parent, Pseq))
	continue;
      
      
      if(d%2!=0)
	numConnectors++;
      numConnectors += descendants(P, T, G, B, Cmax[z], d+1);
      if(d%2==0)
	numConnectors++;
      
#ifdef DEBUG
      int num_debug;
      num_debug = Cmax[z]->seq.size();
      Dist[num_debug]++;
      if(num_debug==1)
	F[Cmax[z]->seq[0]]++;
#endif

    }
    for(int z=0;z<numCmax;z++){
      if(Cmax[z]->I!=NULL)
	delete[] Cmax[z]->I;
      delete Cmax[z];
    }
    Cmax.clear();
  }
  return numConnectors;
}


void initBFSTool(BFSTool B, int n){
  B->marker = 1;
  B->marker_max = MARKER_MAX;
  for(int i=0;i<n;i++){
    B->mark.push_back(0);
    B->bfs_color.push_back(BFS_WHITE);
  }
}
