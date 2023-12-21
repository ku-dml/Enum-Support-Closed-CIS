/********** output stats of instance **********/

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
#include "define.h"
#include "common.h"
#include "mylib.h"
#include "mt19937ar.h"

double cpu_time();

/***** main function *****/
int main(int argc, char *argv[]){
  Param P;
  Tool T;
  Graph G;

  /*** read the arguments ***/
  if(argc<3){
    fprintf(stderr, "usage: %s (pattern_file)(graph_file)\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  P = new struct _Param;
  T = new struct _Tool;
  G = new struct _Graph;
  P->ptn_file = argv[1];
  P->grh_file = argv[2];
  P->reduce = false;

  /*** read the pattern file and the graph file ***/
  preprocess(P, T, G);
  readPattern(P, T, G);
  readGraph(P, T, G);

  cout << G->n << "\t";
  cout << G->items << "\t";

  // edge density
  double er;
  er = 2.0 * G->m / (G->n*(G->n-1));
  cout << G->m << "\t" << er << "\t";

  // item density
  int ir_sum=0,ir_max=-1,ir_min=G->items;
  for(int x=0;x<G->n;x++){
    if(G->V[x]==NULL)
      continue;
    ir_sum += G->V[x]->items;
    if(G->V[x]->items>ir_max)
      ir_max = G->V[x]->items;
    if(G->V[x]->items<ir_min)
      ir_min = G->V[x]->items;
  }
  cout << (double)ir_sum/(double)(G->n*G->items) << "\t";
  cout << (double)ir_min/(double)G->items << "\t";
  cout << (double)ir_max/(double)G->items << "\t";

  // maximum degree
  int max_deg = -1;
  for(int x=0;x<G->n;x++){
    if(G->V[x]==NULL)
      continue;
    int deg = 0;
    AdjList A;
    A = G->V[x]->A;
    while(A!=NULL){
      deg++;
      A = A->next;
    }
    if(deg>max_deg)
      max_deg = deg;
  }
  cout << max_deg << "\n";
  
  /*** postprocess ***/
  delete P;
  delete T;
  delete G;
  
  return EXIT_SUCCESS;
}
