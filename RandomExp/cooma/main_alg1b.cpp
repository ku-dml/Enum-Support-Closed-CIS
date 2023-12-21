/********** A new algorithm for CCIG enumeration problem 
	    (CCIG: Closed Common Itemset Connected subGraph) **********/

/* 
   Implementation notes:

   [usage of indices; there may some exceptions]
   - i,j     ... items
   - x,y,z   ... vertices and edges
   - p,q     ... bitstring
   - b       ... bit
   - k,l,s,t ... indices of others
*/

#define DEBUG
#undef DEBUG

#define STEPWISE
#undef STEPWISE

#define MEM
#undef MEM

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

#include "define.h"
#include "common.h"
#include "mylib.h"
#include "mt19937ar.h"

double cpu_time();

/***** declarations *****/
void algorithm(Param P, Tool T, Graph G, TrieForest &TrieF);
Graph getSubgraphByMark(Graph G, int lower, int upper);
SeqTree computeMaximalFamily(OrderedSeqTree &F, vector<int> &D, int n, int j);


/***** main function *****/
int main(int argc, char *argv[]){
  Param P;
  Tool T;
  Graph G;
  TrieForest TrieF;
  double start_time, fin_time;
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
#ifdef MEM
  printf("BeforeVmSize: %d\nBeforeVmRSS: %d\n", 
	 getVirtualMem(), getPhysicalMem());
#endif 
  
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
  bool safe_terminate = true;
  OrderedSeqTree F,S;
  vector<SeqTree> C_all, C_multi, F_j;
  vector<VertexIDSeq> seq_item;
  vector<int> D, SizeVec;
  VertexIDSeq seq_all, seq, X;
  int x, i, j, p, s, marker, max=-1, size;
  
  double dummy;
  cpu_time();
  dummy = cpu_time();
  
  /*** prepare seq_all ***/
  for(x=0;x<G->n;x++)
    seq_all.push_back(x);

  /*** Step 1: compute C_all and C_multi ***/
  for(i=0;i<G->items;i++){
    // C_all: set of all vertices that have item i
    seq_item.push_back(getSubseqByItem(T, G, i, seq_all));
    marker = markSet(T, G, G, seq_item[i]);
    C_all.push_back(computeSeqTree(T, G, seq_item[i], marker, EVERY_ITEM));
    // update max
    for(auto itr=C_all[i].begin(); itr!=C_all[i].end(); ++itr){
      seq = itr->first;
      size = seq.size();
      if(size > max)
	max = size;
    } 
    // C_multi: subset of C_all such that isolated vertex is not included
    marker = markSet(T, G, G, seq_item[i]);
    C_multi.push_back(computeSeqTree(T, G, seq_item[i], marker, i));
#ifdef DEBUG
    cout << "C_all[" << T->IMapInv[i] <<"] = ";
    outputSeqTree(T, G, C_all[i]);
    cout << "\n";
    cout << "C_multi[" << T->IMapInv[i] << "] =";
    outputSeqTree(T, G, C_multi[i]);
    cout << "\n";
#endif
  }

  // SizeVec: direction of sizewise search
  for(s=max;s>=2;s--)
    SizeVec.push_back(s);

  /*** Step 2 and 3: compute the collections S and F ***/
  for(i=0;i<G->items;i++)
    for(auto itr=C_all[i].begin(); itr!=C_all[i].end(); ++itr){
      seq = itr->first;
      size = seq.size();
      if(size==1)
	S[seq] = itr->second;
      else
	F[seq] = itr->second;
    }
 
  cout << "base_single:\t" << S.size() << "\n";
  cout << "base_multiple:\t" << F.size() << "\n";

  /*** Step 4 to 10: construction of F_j's ***/
#ifdef A1C
  j = 0;
  for(auto itr=F.begin(); itr!=F.end(); ++itr){
    SeqTree tree;
    tree[itr->first] = true;
    F_j.push_back(tree);
#ifdef DEBUG
    cout << "F_" << j << " = ";
    outputSeqTree(T, G, F_j[j]);
    cout << "\n";
#endif
    j++;
  }
#else
  for(x=0;x<G->n;x++)
    D.push_back(-1);
  j = 0;
  while(F.size()>0){
    F_j.push_back(computeMaximalFamily(F, D, G->n, j));
#ifdef DEBUG
    cout << "F_" << j << " = ";
    outputSeqTree(T, G, F_j[j]);
    cout << "\n";
#endif
    for(auto itr=F_j[j].begin(); itr!=F_j[j].end(); ++itr)
      F.erase(itr->first);
    j++;
    if(j==INT_MAX){
      fprintf(stderr, "error: j has become too large.\n");
      exit(EXIT_FAILURE);
    }      
  }
#endif

  p = j;
  
#ifdef DEBUG
  printf("time: %g\n", cpu_time()-dummy);
  printf("p=%d / %d\n", p, G->items);
  dummy = cpu_time();
#endif

  /*** Step 11: initialize TrieF ***/
  for(s=0;s<=max;s++){
    Trie trie;
    trie = new struct _Trie;
    initTrie(trie, s);
    TrieF.push_back(trie);
  }
  for(auto itr=F_j[0].begin(); itr!=F_j[0].end(); ++itr){
    seq = itr->first;
    size = seq.size();
    insertTrie(TrieF[size], seq, Undef, itr->second);
  }
  for(auto itr=S.begin(); itr!=S.end(); ++itr){
    seq = itr->first;
    size = seq.size();
    insertTrie(TrieF[size], seq, Undef, itr->second);
  }
  printf("p:\t%d\n", p);

#ifdef DEBUG
  printf("time: %g\n", cpu_time()-dummy);
#endif

  /*** Step 12-23: for-loop ***/
  for(j=1;j<p;j++){

#ifdef MEM
    if(P->ramub>0 && getPhysicalMem()>P->ramub){
      printf("ABORT_BY_RAMLIMIT\n");
      safe_terminate = false;
      break;
    }
#endif

  
#ifdef DEBUG
    printf("(j=%d/%d)\n", j, p);
    size = 0;
    for(s=0;s<=max;s++){
      cout << " " << TrieF[s]->size;
      size += TrieF[s]->size;
    }
    cout << "  " << size << "\n";
    cout << Skip << "\t" << Swing << "\t" << cpu_time() << "\t" << cpu_time()/(double)size << "\n";
#endif

    // Step 13
    Graph G_j;
    SeqTree::iterator seqitr;
    int lower,upper;
    seqitr = F_j[j].begin();
    lower = markSet(T, G, G, seqitr->first);
    for(; seqitr!=F_j[j].end(); ++seqitr)
      upper = markSet(T, G, G, seqitr->first);
    G_j = getSubgraphByMark(G, lower, upper);
    
    // Step 14: M_temp is generated and M_current is copied to it
    // ---> not needed in the current implementation

    /* Step 15-21: while-loop */
    for(auto s_itr=SizeVec.begin();s_itr!=SizeVec.end();++s_itr){
      s = *s_itr;
      if(s<=P->sigma)
	continue;
      TrieNode u = TrieF[s]->first_leaf;
      while(u!=NULL){
	
	// Step 15 and 16
	if(u->value==j){
#ifdef DEBUG
	  Skip++;
#endif
	  u = u->next;
	  continue;
	}
	X = u->seq;
      
	// Step 17 ---> not needed in the current implementation

	// Step 18
	seq = getIntersection(T, G_j, X);
	marker = markSet(T, G, G_j, seq);
	SeqTree C = computeSeqTree(T, G_j, seq, marker, EVERY_ITEM);
	if(C.empty()){
#ifdef DEBUG
	  Swing++;
#endif
	  u = u->next;
	  continue;   
	} 
      
	// Step 19
	for(auto itr=C.begin(); itr != C.end(); ++itr){
	  seq = itr->first;
	  size = seq.size();
	  if(size>=P->sigma)
	    insertTrie(TrieF[size], seq, j, itr->second);
	}
	// Step 20 ---> not needed in the current implementation

	u = u->next;
      }
    }
    
    // Step 22
    for(auto itr=F_j[j].begin(); itr!=F_j[j].end(); ++itr){
      seq = itr->first;
      size = seq.size();
      if(size>=P->sigma)
	insertTrie(TrieF[size], seq, j, itr->second);
    }
    
    deleteSubgraph(G_j);
    delete G_j;

#ifdef STEPWISE
    size = 0;
    for(s=0;s<=max;s++)
      size += TrieF[s]->size;
    cout << "step\t" << j << "\t" << size << "\t" << cpu_time()-dummy << "\n";
#endif
    
    if(P->tlim>0. && cpu_time()>P->tlim){
      printf("ABORT_BY_TIMELIMIT\n");
      safe_terminate = false;
      break;
    }
  }
  if(safe_terminate)
    printf("TERMINATE_SAFELY\n");
  
#ifdef DEBUG
  size = 0;
  for(s=0;s<=max;s++){
    cout << " " << TrieF[s]->size;
    size += TrieF[s]->size;
  }
  cout << " " << size << "\n";
  cout << Skip << "\t" << Swing << "\t" << cpu_time() << "\n";
#endif

  
  if(P->theta<=1 && P->delta<0.000001)
    return;
    
  for(s=0;s<=max;s++){
    TrieF[s]->size = 0;
    if(TrieF[s]->h < P->sigma){
      TrieF[s]->first_leaf = NULL;
      continue;
    }
    TrieNode u = TrieF[s]->first_leaf;
    TrieNode prev = NULL;
    while(u!=NULL){
      if(isIntersectionSmall(T->B, T->d, G, u->seq, P->theta) ||
	 getDensity(u->seq.size(), u->deg/2) < P->delta){
	if(prev==NULL)
	  TrieF[s]->first_leaf = u->next;
	else
	  prev->next = u->next;
      }
      else{
	prev = u;
	TrieF[s]->size++;
      }
      u = u->next;
    }
  }
  
#if 0
  for(s=0;s<=max;s++){
    if(TrieF[s]->h < P->sigma)
      continue;
    TrieNode u = TrieF[s]->first_leaf;
    while(u!=NULL){
      //outputSet(T, G, u->seq);
      printf("\t%d %d %d\t%g\n", u->seq.size(), getNumItems(T->B, T->d, G, u->seq), u->deg, getDensity(u->seq.size(), u->deg/2));
      u = u->next;
    }
  }
#endif
}



/***** get the subgraph of G induced by vertices
       whose marks are in [lower,upper];
       edges are drawn only for two vertices that have the same mark *****/  
Graph getSubgraphByMark(Graph G, int lower, int upper){
  Graph sub;
  int x,k,size;
  // construction of vertices
  sub = new struct _Graph;
  sub->V = new Vertex[G->n];
  for(x=0;x<G->n;x++){
    if(G->V[x]==NULL){
      sub->V[x] = NULL;
      continue;
    }
    if(G->V[x]->mark<lower || G->V[x]->mark>upper){
      sub->V[x] = NULL;
      continue;
    }
    sub->V[x] = new struct _Vertex;
    sub->V[x]->id = G->V[x]->id;
    sub->V[x]->I = G->V[x]->I;
    sub->V[x]->A = NULL;
    sub->V[x]->mark = G->V[x]->mark;
    sub->V[x]->bfs_color = G->V[x]->bfs_color;
    sub->V[x]->bfs_conn = G->V[x]->bfs_conn;
  }
  sub->E = NULL; // pending (may not be used for the time being)
  sub->n = G->n;
  sub->m = G->m; // pending (may not be used for the time being)
  sub->items = G->items;

  // construction of edges
  for(x=0;x<G->n;x++){
    if(sub->V[x] == NULL)
      continue;
    AdjList l;
    l = G->V[x]->A;
    while(l!=NULL){
      int z = l->v->id;
      if(sub->V[z]!=NULL && sub->V[x]->mark==sub->V[z]->mark){
	AdjList newl;
	newl = new struct _AdjList;
	newl->v = sub->V[z];
	newl->e = new struct _Edge;
	newl->e->V[0] = sub->V[x];
	newl->e->V[1] = sub->V[z];
	newl->e->w = l->e->w;
	if(sub->V[x]->A==NULL){
	  sub->V[x]->A = newl;
	  newl->next = NULL;
	}
	else{
	  newl->next = sub->V[x]->A;
	  sub->V[x]->A = newl;
	}
      }
      l = l->next;
    }
  }
  // postprocess
  return sub;
}



/***** compute a maximal family of disjoint sets from F *****/
SeqTree computeMaximalFamily(OrderedSeqTree &F, vector<int> &D, int n, int j){
  SeqTree F_j;
  VertexIDSeq seq;
#ifdef COLMAX
  OrderedSeqTree::reverse_iterator itr;
  // largest subset is inserted
  itr = F.rbegin();
  seq = itr->first;
  F_j[seq] = itr->second;
  for(auto i=seq.begin(); i!=seq.end(); ++i)
    D[*i] = j;

  // search for subsets to be inserted
  for(; itr!=F.rend(); ++itr){
    bool flag = true;
    seq = itr->first;
    for(auto i=seq.begin(); i!=seq.end(); ++i)
      if(D[*i] == j){
	flag = false;
	break;
      }
    if(!flag)
      continue;
    F_j[seq] = itr->second;
    for(auto i=seq.begin(); i!=seq.end(); ++i)
      D[*i] = j;    
  }
#else
  OrderedSeqTree::iterator itr;
  // smallest subset is inserted
  itr = F.begin();
  seq = itr->first;
  F_j[seq] = itr->second;
  for(auto i=seq.begin(); i!=seq.end(); ++i)
    D[*i] = j;

  // search for subsets to be inserted
  for(; itr!=F.end(); ++itr){
    bool flag = true;
    seq = itr->first;
    for(auto i=seq.begin(); i!=seq.end(); ++i)
      if(D[*i] == j){
	flag = false;
	break;
      }
    if(!flag)
      continue;
    F_j[seq] = itr->second;
    for(auto i=seq.begin(); i!=seq.end(); ++i)
      D[*i] = j;    
  }
#endif
  return F_j;
}


