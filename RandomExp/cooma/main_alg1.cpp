/********** A new algorithm for CCIG enumeration problem 
	    (CCIG: Closed Common Itemset Connected subGraph) **********/

/* 
   Implementation notes:

   [usage of indices; there may be some exceptions]
   - i,j     ... items
   - x,y,z   ... vertices and edges
   - p,q     ... bitstring
   - b       ... bit
   - k,l,s,t ... indices of others
*/

#define DEBUG
#undef DEBUG

#define STEPWISE

#define MEM
//#undef MEM

/***** includes *****/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

// C++ libraries
#include <iostream>
#include <unordered_map> // hash
#include <map>           // red-black tree
#include <queue>
#include <vector>
#include <algorithm>
using namespace std;

// custom libraries
#ifdef MEM
#include "mem.cpp"
#endif

// custom libraries
#include "define.h"
#include "common.h"
#include "mylib.h"
#include "mt19937ar.h"

double cpu_time();

/***** declarations *****/
void algorithm(Param P, Tool T, Graph G, TrieForest &TrieF);
int chooseInitItem(Param P, Graph G);

/***** main function *****/
int main(int argc, char *argv[]){
  Param P;
  Tool T;
  Graph G;
  TrieForest TrieF;
  double start_time, fin_time;
  int s;
  unsigned long int size = 0;
  unsigned long int vertex_sum = 0;
  unsigned long int edge_sum = 0;

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

  printf("vertices:\t%d\n", G->n);
  printf("edges:\t%d\n", G->m);
  printf("edge_density:\t%g\n", getDensity(G->n, G->m));
  printf("items:\t%d\n", G->items);
  printf("item_density:\t%g\n", G->item_density);
  
#ifdef DEBUG
  outputGraph(T, G);
#endif

  /*** run algorithm ***/
  start_time = cpu_time();
  algorithm(P, T, G, TrieF);
  fin_time = cpu_time();
  
#ifdef MEM
  printf("VIRTUAL: %d\nPHYSICAL: %d\n", 
	 getVirtualMem(), getPhysicalMem());
#endif 
  
  for(auto itr=TrieF.begin(); itr!=TrieF.end(); itr++){
    if((*itr)->h < P->sigma)
      continue;
    TrieNode u = (*itr)->first_leaf;
    while(u!=NULL){
      size++;
      vertex_sum += u->seq.size();
      edge_sum += u->deg/2;

      /*
      outputSet(T,G,u->seq);
      cout << " " << u->deg/2 << "\n";
      */
      
      u = u->next;
    }
  }
    
#ifdef DEBUG
  cout << "M_current:\n";
  outputTrieForest(P, T, G, TrieF);
#endif
  cout << "subgraphs:\t" << size << "\n";
  cout << "vertex_sum:\t" << vertex_sum << "\n";
  cout << "edge_sum:\t" << edge_sum << "\n";
  cout << "cpu_time:\t" << fin_time-start_time << "\n";
  if(P->reduce){
    printf("reduced_vertices:\t%d\n", T->reduced_vertices);
    printf("reduced_edges:\t%d\n", T->reduced_edges);
  }

  //outputTrieForest(P, T, G, TrieF);

  if(P->distname!=NULL)
    outputDist(fopen(P->distname, "w"), T, G, TrieF);
    
  /*** postprocess ***/
  delete P;
  delete T;
  delete G;
  
  return EXIT_SUCCESS;
}

#ifdef DEBUG
int Skip=0;  // # of skipped components
int Swing=0; // # of X's s.t. X and G_i are disjoint
#endif

/***** main algorithm *****/
void algorithm(Param P, Tool T, Graph G, TrieForest &TrieF){
  vector<SeqTree> C, F; // elementary vertex sets for each item 
  vector<VertexIDSeq> seq_item; // all vertices that have item i
  VertexIDSeq seq_all, seq, X;
  vector<int> SizeVec;
  int i, i_1, s, marker, max=-1, size;

  double dummy;
  cpu_time();
  dummy = cpu_time();
  
  /*** prepare seq_all ***/
  for(i=0;i<G->n;i++)
    seq_all.push_back(i);

  /*** Step 1: compute C and F ***/
  for(i=0;i<G->items;i++){
    // C: set of all vertices that have item i
    seq_item.push_back(getSubseqByItem(T, G, i, seq_all));
    marker = markSet(T, G, G, seq_item[i]);
    C.push_back(computeSeqTree(T, G, seq_item[i], marker, EVERY_ITEM));
    // update max: # of vertices in the largest component in C
    for(auto itr=C[i].begin(); itr!=C[i].end(); ++itr){
      seq = itr->first;
      size = seq.size();
      if(size > max)
	max = size;
    }
    // F: subset of C such that isolated vertices are not included
    marker = markSet(T, G, G, seq_item[i]);
    F.push_back(computeSeqTree(T, G, seq_item[i], marker, i));
#ifdef DEBUG
    cout << "C[" << T->IMapInv[i] <<"] = ";
    outputSeqTree(T, G, C[i]);
    cout << "\n";
    cout << "F[" << T->IMapInv[i] << "] =";
    outputSeqTree(T, G, F[i]);
    cout << "\n";
#endif
  }

  // SizeVec: direction of sizewise search
  for(s=max;s>=2;s--)
    SizeVec.push_back(s);
  
  /*** Step 2: choose i_1 ***/
  i_1 = chooseInitItem(P, G);
#ifdef DEBUG
  printf("\nItem %d is chosen as i_1\n\n", T->IMapInv[i_1]);
#endif
  
  /*** Step 3: initialize TrieF ***/
  for(s=0;s<=max;s++){
    Trie trie;
    trie = new struct _Trie;
    initTrie(trie, s);
    TrieF.push_back(trie);
  }
  for(i=0;i<G->items;i++){
    for(auto itr=C[i].begin(); itr != C[i].end(); ++itr){
      seq = itr->first;
      size = seq.size();
      if(i==i_1 || size==1)
	insertTrie(TrieF[size], seq, Undef, itr->second);
    }
  }
  
  /*** Step 4-15: for-loop ***/
  for(i=0;i<G->items;i++){
  
#ifdef DEBUG
    printf("(i=%d <%d>/%d)\n", i, T->IMapInv[i], G->items);
    size = 0;
    for(s=0;s<=max;s++){
      cout << " " << TrieF[s]->size;
      size += TrieF[s]->size;
    }
    cout << "  " << size << "\n";
    cout << "Skip="<< Skip << "\tSwing=" << Swing << "\tcpu_time=" << cpu_time() << "\ttime_per_set=" << cpu_time()/(double)size<< "\n";
#endif
    
    if(i==i_1)
      continue;

    // Step 5
    Graph G_i = getSubgraph(G, seq_item[i]);

    /* Step 7-13: while-loop */
    for(auto s_itr=SizeVec.begin();s_itr!=SizeVec.end();++s_itr){
      s = *s_itr;
      TrieNode u = TrieF[s]->first_leaf;
      while(u!=NULL){
	/***
	    NOTE:
	    - u->seq is X chosen in Step 8.
	    - u->value is the latest iteration when the key is searched
	***/
	if(u->value==i){
#ifdef DEBUG
	  Skip++;
#endif
	  u = u->next;
	  continue;
	}

	// Step 8: set treeitr->first to X
	X = u->seq;

	// Step 9: you do not need to take it into account in this implementation

	// Step 10
	seq = getIntersection(T, G_i, X);
	marker = markSet(T, G, G_i, seq);
	SeqTree C = computeSeqTree(T, G_i, seq, marker, i);
	
	if(C.empty()){
#ifdef DEBUG
	  Swing++;
#endif
	  u = u->next;
	  continue;    
	}	

	// Step 11
	for(auto itr=C.begin(); itr != C.end(); ++itr){
	  seq = itr->first;
	  size = seq.size();
	  insertTrie(TrieF[size], seq, i, itr->second);
	}
     	// Step 12: you do not need to take it into account in this implementation

	u = u->next;
      }
    }
    // Step 14
    for(auto itr=F[i].begin(); itr != F[i].end(); ++itr){
      seq = itr->first;
      size = seq.size();
      insertTrie(TrieF[size], seq, i, itr->second);
    }
    deleteSubgraph(G_i);
    delete G_i;
    
#ifdef STEPWISE
    size = 0;
    for(s=0;s<=max;s++)
      size += TrieF[s]->size;
    cout << "step\t" << i << "\t" << size << "\t" << cpu_time()-dummy << "\n";
#endif
    
    if(P->tlim>0. && cpu_time()>P->tlim){
      printf("ABORT_BY_TIMELIMIT\n");
      break;
    }
  }
#ifdef DEBUG
  size = 0;
  for(s=0;s<=max;s++){
    cout << " " << TrieF[s]->size;
    size += TrieF[s]->size;
  }
  cout << " " << size << "\n";
  cout << Skip << "\t" << Swing << "\t" << cpu_time() << "\n";
#endif
  
  if(P->theta>1){
    for(s=0;s<=max;s++){
      TrieF[s]->size = 0;
      TrieNode u = TrieF[s]->first_leaf;
      while(u!=NULL){
	if(!isIntersectionSmall(T->B, T->d, G, u->seq, P->theta))
	  TrieF[s]->size++;
	u = u->next;
      }
    }
  }
}


/***** Step2 of the algorithm: choose i_1 *****/
int chooseInitItem(Param P, Graph G){
  int i_1;
  i_1 = (int)(genrand_real2()*(double)G->items);
  return i_1;
}



