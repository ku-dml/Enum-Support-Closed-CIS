/********** common.cpp **********/

/***** includes *****/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>

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


/***** check the arguments *****/
void checkArgs(int argc, char *argv[]){
  if(argc<4){
    fprintf(stderr, "usage: %s (pattern_file)(graph_file)(theta)[(param1)(value1)(param2)...]\n\n", argv[0]);
    fprintf(stderr, "list of params:\n");
    fprintf(stderr, "  -seed <INT> ... random seed (%d)\n", INI_seed);
    fprintf(stderr, "  -sigma <INT> ... minimum size of connector (%d)\n", INI_sigma);
    fprintf(stderr, "  -delta <DOUBLE> ... minimum density of a connector (%g)\n", INI_delta);
    fprintf(stderr, "  -tlim <DOUBLE> ... time limit in sec; if <0, it is inf. (%g)\n", INI_tlim);
    fprintf(stderr, "  -reduce <BOOL> ... whether reduction is invoked (%d)\n", INI_reduce);
    fprintf(stderr, "  -vtable <STR> ... file that describes relation between vertex id & name; itable should be also specified\n");
    fprintf(stderr, "  -itable <STR> ... file that describes relation between item id & name; vtable should be also specified\n");
    fprintf(stderr, "  -outname <STR> ... filename of output from vtable & itable (%s)\n", INI_outname);
    fprintf(stderr, "  -distname <STR> ... filename for distribution\n", INI_outname);
    fprintf(stderr, "  -ramub <INT> ... upper bound on used amount of RAM (%d)\n", INI_ramub);
    exit(EXIT_FAILURE);
  }
}


/***** read the arguments *****/
void readArgs(Param P, Graph G, int argc, char *argv[]){
  char *arg;
  int k;
  P->ptn_file = argv[1];
  P->grh_file = argv[2];
  P->theta = atoi(argv[3]);
  P->seed = INI_seed;
  P->sigma = INI_sigma;
  P->delta = INI_delta;
  P->tlim = DBL_MAX; // negative means infinity
  P->reduce = INI_reduce;
  P->vtable = NULL;
  P->itable = NULL;
  P->outname = new char[LEN_FILENAME];
  P->distname = NULL;
  P->ramub = INI_ramub;
  strcpy(P->outname, INI_outname);
  
  for(k=4;k<argc;k+=2){
    if(k==argc-1 || argv[k][0] != '-'){
      fprintf(stderr, "error: arguments are not set appropriately.\n");
      exit(EXIT_FAILURE);
    }
    arg = argv[k]+1;
    if(strcmp(arg, "seed")==Equiv)
      P->seed = atoi(argv[k+1]);
    else if(strcmp(arg, "sigma")==Equiv)
      P->sigma = atoi(argv[k+1]);
    else if(strcmp(arg, "delta")==Equiv)
      P->delta = atof(argv[k+1]);
    else if(strcmp(arg, "tlim")==Equiv)
      P->tlim = atof(argv[k+1]);
    else if(strcmp(arg, "reduce")==Equiv)
      P->reduce = (bool)atoi(argv[k+1]);
    else if(strcmp(arg, "vtable")==Equiv)
      P->vtable = argv[k+1];
    else if(strcmp(arg, "itable")==Equiv)
      P->itable = argv[k+1];
    else if(strcmp(arg, "outname")==Equiv)
      strcpy(P->outname, argv[k+1]);
    else if(strcmp(arg, "distname")==Equiv){
      P->distname = new char[LEN_FILENAME];
      strcpy(P->distname, argv[k+1]);
    }
    else if(strcmp(arg, "ramub")==Equiv)
      P->ramub = atoi(argv[k+1]);
    else{
      fprintf(stderr, "error: <%d> is illegal parameter.\n", arg);
      exit(EXIT_FAILURE);
    }
  }
  if((P->vtable==NULL && P->itable!=NULL) ||
     (P->vtable!=NULL && P->itable==NULL)){
    fprintf(stderr, "error: if one of vtable and itable is specified, then the other should be also specified.\n");
    exit(EXIT_FAILURE);
  }  
}


/***** preprocess *****/
void preprocess(Param P, Tool T, Graph G){
  char *str,*str_id,*str_item,**s;
  FILE *in;
  double w;
  int b,i,j,x,y,items;

  str = new char[LINE_MAX]; 
  str_id = new char[LINE_MAX];
  str_item = new char[LINE_MAX];

  /* first scan of pattern file */
  G->n = 0;
  G->m = 0;
  G->items = 0;
  G->item_density = 0.;
  in = open_file(P->ptn_file, "r");
  while(fgets(str, LINE_MAX-1, in)!=NULL){
    if(str[0]=='#')
      continue;
    sscanf(str, "%s%s", str_id, str_item);
    auto itr = T->VMap.find(atoi(str_id));
    if(itr==T->VMap.end()){ // if str_id is not in the VMap
      T->VMap[atoi(str_id)] = G->n;
      T->VMapInv[G->n] = atoi(str_id);
      G->n++;
    }
    s = split(str_item, ',');
    items = getNumChar(str_item, ',')+1;
    for(j=0;j<items;j++){
      i = atoi(s[j]);
      auto itr = T->IMap.find(i);
      if(itr==T->IMap.end()){ // if i is not in the IMap
	T->IMap[i] = G->items;
	T->IMapInv[G->items] = i;
	G->items++;
      }
      free(s[j]);
    }
    free(s);
  }
  fclose(in);

  // G->V
  G->V = new Vertex[G->n];
  
  /* first scan of graph file */
  in = open_file(P->grh_file, "r");
  while(fgets(str, LINE_MAX-1, in)!=NULL){
    if(str[0]=='#')
      continue;
    sscanf(str, "%d%lf%d", &x, &w, &y);
    if(x!=y)
      G->m++;
  }

  // G->E
  G->E = new Edge[G->m];

  /* initialize T */
  T->B[0] = 1;
  for(b=1;b<BIT_LENGTH;b++)
    T->B[b] = T->B[b-1] << 1;
  T->d = G->items/BIT_LENGTH + 1;
  T->marker = 1;
  T->marker_max = MARKER_MAX;
  T->reduced_vertices = 0;
  T->reduced_edges = 0;

  delete[] str;
  delete[] str_id;
  delete[] str_item;
}


/***** read the pattern file *****/
void readPattern(Param P, Tool T, Graph G){
  char *str,*str_id,*str_item,**s;
  FILE *in;
  double item_density = 0.;
  int x,i,j,p;

  str = new char[LINE_MAX]; 
  str_id = new char[LINE_MAX];
  str_item = new char[LINE_MAX];
  in = open_file(P->ptn_file, "r");

  while(fgets(str,LINE_MAX-1,in)!=NULL){
    if(str[0]=='#')
      continue;
    strcpy(str_item, ""); // for the case itemset is empty
    sscanf(str, "%s%s", str_id, str_item);
    x = T->VMap[atoi(str_id)];
    G->V[x] = new struct _Vertex;
    G->V[x]->id = x;
    G->V[x]->I = new BitString[T->d];
    for(p=0;p<T->d;p++)
      G->V[x]->I[p] = 0;
    G->V[x]->A = NULL;
    G->V[x]->mark = 0;
    G->V[x]->bfs_color = BFS_WHITE;
    if(strcmp(str_item, "")==Equiv || strcmp(str_item, "!")==Equiv)
      G->V[x]->items = 0;
    else{
      s = split(str_item, ',');
      G->V[x]->items = getNumChar(str_item, ',')+1;
      item_density += (double)(G->V[x]->items);
      for(j=0;j<G->V[x]->items;j++){
	i = T->IMap[atoi(s[j])];
	addItem(T->B, T->d, G->V[x]->I, i);
	free(s[j]);
      }
      free(s);
    }
    /* Reduction 1 */
    if(P->reduce && G->V[x]->items<P->theta){
      delete[] G->V[x]->I;
      delete G->V[x];
      G->V[x] = NULL;
      T->reduced_vertices++;
    }
  }

  G->item_density = (item_density/(double)G->items) / (double)G->n;
  
  delete[] str;
  delete[] str_id;
  delete[] str_item;
  fclose(in);
}


/***** read the graph file *****/
void readGraph(Param P, Tool T, Graph G){
  char *str;
  FILE *in;
  double w;
  int a,b,x,y,z=0,line=0;
  
  str = new char[LINE_MAX];
  in = open_file(P->grh_file, "r");

  while(fgets(str,LINE_MAX-1,in)!=NULL){
    line++;
    if(str[0]=='#') // comment is skipped
      continue;
    sscanf(str, "%d%lf%d", &a, &w, &b);
    if(a==b){
      fprintf(stderr, "The line %d is ignored since %d=%d (self-loop).\n", line, a, b);
      continue;
    }
    x = T->VMap[a];
    y = T->VMap[b];
    if(P->reduce){
      if(G->V[x]==NULL || G->V[y]==NULL)
	continue;
      /* reduction 2 */
      ItemStr I;
      int p,items=0;
      I = new BitString[T->d];
      for(p=0;p<T->d;p++)
	I[p] = G->V[x]->I[p] & G->V[y]->I[p];
      for(p=0;p<T->d;p++)
	items += WP3(I[p]);
      delete[] I;
      if(items<P->theta){
	T->reduced_edges++;
	continue;
      }
    }
    G->E[z] = new struct _Edge;
    G->E[z]->V[0] = G->V[x];
    G->E[z]->V[1] = G->V[y];
    G->E[z]->w = w;

    AdjList list;
    bool added = true;

    list = new struct _AdjList;
    list->v = G->E[z]->V[1];
    list->e = G->E[z];

    //#define VERSION_20180225
#ifdef VERSION_20180225
    if(G->E[z]->V[0]->A==NULL){
      G->E[z]->V[0]->A = list;
      list->next = NULL;
    }
    else{
      list->next = G->E[z]->V[0]->A;
      G->E[z]->V[0]->A = list;
    }
#else
    if(G->V[x]->A==NULL){
      G->V[x]->A = list;
      list->next = NULL;
    }
    else{
      if(list->v->id > G->V[x]->A->v->id){
	list->next = G->V[x]->A;
	G->V[x]->A = list;
      }
      else{
	AdjList prev,l;
	bool twice=false;
	prev = G->V[x]->A;
	l = prev->next;
	while(l!=NULL){
	  if(list->v->id == l->v->id){
	    fprintf(stderr, "warning: edge (%d,%d) is found twice.\n", G->V[x]->id, l->v->id);
	    twice = true;
	    break;
	  }
	  else if(list->v->id > l->v->id)
	    break;
	  prev = l;
	  l = l->next;
	}
	if(twice==true){
	  delete G->E[z];
	  delete list;
	  added = false;
	}
	else{
	  list->next = prev->next;
	  prev->next = list;
	}
      }
    }
#endif
    if(added)
      z++;
  }

  G->m = countNumEdges(G);

  delete[] str;
  fclose(in);
}

/***** count the number of edges *****/
int countNumEdges(Graph G){
  int m=0;
  for(int x=0;x<G->n;x++){
    if(G->V[x]==NULL)
      continue;
    int deg=0;
    for(AdjList list=G->V[x]->A;list!=NULL;list=list->next)
      deg++;
    m += deg;
  }
  if(m%2==1)
    fprintf(stderr, "warning: degree sum is not even.\n");
  return m/2;
}

/***** add item i to I, which consists of d BitString *****/
void addItem(BitString B[BIT_LENGTH], int d, ItemStr I, int i){
  int p;
  p = i/BIT_LENGTH;
  I[p] += B[i%BIT_LENGTH];
}


/***** return whether I has item i *****/
bool hasItem(BitString B[BIT_LENGTH], int d, ItemStr I, int i){
  int p;
  p = i/BIT_LENGTH;
  if(B[i%BIT_LENGTH]&I[p])
    return true;
  return false;
}


/***** return whether the size of the union of itemsets over seq is < theta *****/
bool isUnionSmall(BitString B[BIT_LENGTH], int d, Graph G,
		  VertexIDSeq seq, int theta){
  bool flag = true;
  ItemStr I;
  int p,x,size;
  I = new BitString[d];
  for(p=0;p<d;p++) I[p] = 0;
  size = seq.size();
  for(x=0;x<size;x++){
    if(G->V[seq[x]]->items >= theta){
      flag = false;
      break;
    }
    for(p=0;p<d;p++)
      I[p] = I[p] | G->V[seq[x]]->I[p];
  }
  if(flag){
    int items = 0;
    for(p=0;p<d;p++)
      items += WP3(I[p]);
    if(items>=theta)
      flag = false;
  }
  delete[] I;
  return flag;
}


/***** return whether the size of the intersection of itemsets over seq is < theta *****/
bool isIntersectionSmall(BitString B[BIT_LENGTH], int d, Graph G,
			 VertexIDSeq seq, int theta){
  bool flag = false;
  ItemStr I;
  int p,x,size;
  I = new BitString[d];
  for(p=0;p<d;p++) I[p] = (~(BitString)0);
  size = seq.size();
  for(x=0;x<size;x++){
    if(G->V[seq[x]]->items < theta){
      flag = true;
      break;
    }
    for(p=0;p<d;p++)
      I[p] = I[p] & G->V[seq[x]]->I[p];
  }
  if(!flag){
    int items = 0;
    for(p=0;p<d;p++)
      items += WP3(I[p]);
    if(items<theta)
      flag = true;
  }
  delete[] I;
  return flag;
}

int getNumItems(BitString B[BIT_LENGTH], int d, Graph G, VertexIDSeq seq){
  ItemStr I;
  int p,x,size,items=0;
  I = new BitString[d];
  for(p=0;p<d;p++) I[p] = (~(BitString)0);
  size = seq.size();
  for(x=0;x<size;x++){
    for(p=0;p<d;p++)
      I[p] = I[p] & G->V[seq[x]]->I[p];
  }
  for(p=0;p<d;p++)
    items += WP3(I[p]);
  delete[] I;
  return items;
}


ItemIDSeq getItemIntersection(Tool T, Graph G, VertexIDSeq seq){
  VertexIDSeq itmseq;
  ItemStr I;
  int i,p,x;
  I = new BitString[T->d];
  for(p=0;p<T->d;p++) I[p] = (~(BitString)0);
  for(auto itr=seq.begin(); itr!=seq.end(); ++itr){
    x = *itr;
    for(p=0;p<T->d;p++)
      I[p] = I[p] & G->V[x]->I[p];
  }
  for(i=0;i<G->items;i++)
    if(hasItem(T->B, T->d, I, i))
      itmseq.push_back(i);
  return itmseq;
}


/***** get subsequence of seq such that 
       the vertices have item i *****/
VertexIDSeq getSubseqByItem(Tool T, Graph G, int i, VertexIDSeq seq){
  VertexIDSeq subseq;
  int x,k,size;
  size=seq.size();
  for(k=0;k<size;k++){
    x = seq[k];
    if(G->V[x]!=NULL && hasItem(T->B, T->d, G->V[x]->I, i))
      subseq.push_back(x);
  }
  return subseq;
}


/***** get the subseq of seq such that subseq is in G *****/
VertexIDSeq getIntersection(Tool T, Graph G, VertexIDSeq seq){
  VertexIDSeq subseq;
  int x;  
  int size=seq.size(),k;
  for(k=0;k<size;k++){
    x = seq[k];
    if(G->V[x]!=NULL)
      subseq.push_back(x);
  }
  return subseq;
}



/***** get the subgraph of G induced by vertices in seq *****/  
Graph getSubgraph(Graph G, VertexIDSeq seq){
  Graph sub;
  int x,k,size;
  // preprocess
  sub = new struct _Graph; //(Graph)mallocE(sizeof(struct _Graph));
  sub->V = new Vertex[G->n]; //(Vertex*)mallocE(G->n*sizeof(Vertex));
  for(x=0;x<G->n;x++)
    sub->V[x] = NULL;
  sub->E = NULL; // pending (may not be used for the time being)
  sub->n = G->n;
  sub->m = G->m; // pending (may not be used for the time being)
  sub->items = G->items;

  // construction of vertices
  size=seq.size();
  for(k=0;k<size;k++){
    x = seq[k];
    if(G->V[x]==NULL){
      fprintf(stderr, "error: invalid subsetting the graph\n");
      exit(EXIT_FAILURE);
    }
    sub->V[x] = new struct _Vertex;
    sub->V[x]->id = G->V[x]->id;
    sub->V[x]->I = G->V[x]->I;
    sub->V[x]->A = NULL;
    sub->V[x]->mark = G->V[x]->mark;
    sub->V[x]->bfs_color = G->V[x]->bfs_color;
    sub->V[x]->bfs_conn = G->V[x]->bfs_conn;
  }

  // construction of edges
  for(k=0;k<size;k++){
    x = seq[k];
    AdjList l;
    l = G->V[x]->A;
    while(l!=NULL){
      int z = l->v->id;
      if(sub->V[z]!=NULL){
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


/***** delete a Graph object *****/
void deleteSubgraph(Graph G){
  int x;
  for(x=0;x<G->n;x++){
    if(G->V[x]==NULL)
      continue;
    AdjList l = G->V[x]->A;
    while(l!=NULL){
      AdjList del = l;
      l = l->next;
      delete del->e;
      delete del;
    }
    delete G->V[x];
  }
  delete[] G->V;
}


/***** mark the vertices in the seq; returns the used marker *****/
int markSet(Tool T, Graph U, Graph G, VertexIDSeq seq){
  int k,size,marker;
  marker = T->marker;
  increaseMarker(T, U, G);
  size=seq.size();
  for(k=0;k<size;k++)
    G->V[seq[k]]->mark = marker;
  return marker;
}


/***** increase marker by 1; returns whether reset is done *****/
bool increaseMarker(Tool T, Graph U, Graph G){
  bool flag=false;
  int x;
  T->marker++;
  // reset the marks of all vertices if marker=marker_max
  if(T->marker>T->marker_max){
    for(x=0;x<U->n;x++)
      U->V[x]->mark = 0;
    T->marker = 1;
    flag = true;
  }
  return flag;
}


/***** BFS that starts from vertex v in S
       conn = id of connected component *****/
int markedBFS(Tool T, Graph G, int marker, Vertex v, int conn){
  AdjList A;
  Vertex u;
  queue<Vertex> Q;
  //int mindeg = G->n;
  int degsum = 0;

  Q.push(v);
  v->bfs_color = BFS_GRAY;
  while(Q.empty()==false){
    int deg=0;
    u = Q.front();
    u->bfs_color = BFS_BLACK;
    u->bfs_conn = conn;
    Q.pop();
    A = u->A;
    while(A!=NULL){
      if(A->v->mark == marker){
	deg++;
	if(A->v->bfs_color==BFS_WHITE){
	  Q.push(A->v);
	  A->v->bfs_color = BFS_GRAY;
	}
      }
      A = A->next;
    }
    /*
    if(deg<mindeg)
      mindeg = deg;
    */
    degsum += deg;
  }
  //return mindeg;
  return degsum;
  
}


/***** compute connected components of (S,E(i));
       i is either item index or macro EVERY_ITEM.
       In the latter case, all edges are considered *****/
SeqTree computeSeqTree(Tool T, Graph G, VertexIDSeq seq, int marker, int i){
  SeqTree SeqT;
  AdjList A;
  vector<VertexIDSeq> S;
  vector<int> degs;
  int x,k,t,size,conn=0;
  // initialize color
  size=seq.size();
  for(k=0;k<size;k++){
    x = seq[k];  
    G->V[x]->bfs_conn = Undef;
    if(i==EVERY_ITEM)
      G->V[x]->bfs_color = BFS_WHITE;
    else{ // if i!=EVERY_ITEM, x is considered only if it has a neighbor that has item i
      G->V[x]->bfs_color = BFS_BLACK;
      A = G->V[x]->A;
      while(A!=NULL){
	if(hasItem(T->B, T->d, A->v->I, i)){
	  G->V[x]->bfs_color = BFS_WHITE;
	  break;
	}
	A = A->next;
      }
    }
  }
  // successive BFS
  for(k=0;k<size;k++){
    x = seq[k];
    if(G->V[x]->bfs_color == BFS_WHITE){
      int deg;
      deg = markedBFS(T, G, marker, G->V[x], conn);
      degs.push_back(deg);
      conn++;
    }
  }
  // construct SeqT
  if(conn==0){
    SeqT.clear();
    return SeqT;
  }
  for(k=0;k<conn;k++){
    VertexIDSeq s;
    S.push_back(s);
  }
  for(t=0;t<size;t++){
    x = seq[t];
    k = G->V[x]->bfs_conn;
    if(k!=Undef)
      S[k].push_back(G->V[x]->id);
  }
  for(k=0;k<conn;k++)
    SeqT[S[k]] = degs[k];

  return SeqT;
}


/***** output seq to stdout *****/
void outputSet(Tool T, Graph G, VertexIDSeq seq){
  int x;
  cout << "{";
  for(auto itr=seq.begin(); itr!=seq.end(); ++itr){
    x = *itr;
    cout << " " << T->VMapInv[G->V[x]->id];
  }
  cout << "}";
}


void outputSeqTree(Tool T, Graph G, SeqTree &SeqT){
  for(auto itr=SeqT.begin(); itr!=SeqT.end(); ++itr)
    outputSet(T, G, itr->first);
}

 
/***** output SeqF to stdout *****/
void outputSeqForest(Tool T, Graph G, SeqForest &SeqF){
  for(auto treeitr=SeqF.rbegin(); treeitr!=SeqF.rend(); ++treeitr){
    SeqTree tree = *treeitr;
    int size;
    size = tree.size();
    if(size > 0)
      outputSeqTree(T, G, tree);
  }
  cout << "\n";
}


/***** output solution in the standard form *****/
void outputInStdform(Tool T, Graph G, SeqForest &SeqF){
  OrderedSeqTree RB;
  VertexIDSeq seq;
  int x; 
  for(auto treeitr=SeqF.begin(); treeitr!=SeqF.end(); ++treeitr){
    SeqTree tree = *treeitr;
    for(auto itr=tree.begin(); itr!=tree.end(); ++itr){
      seq = itr->first;
      sort(seq.begin(), seq.end());
      RB[seq] = true;
    }
  }
  for(auto itr=RB.begin(); itr!=RB.end(); ++itr){
    seq = itr->first;
    cout << seq.size();
    for(auto itr2=seq.begin(); itr2!=seq.end(); ++itr2){
      x = *itr2;
      cout << " " << T->VMapInv[G->V[x]->id];
    }
    cout << "\n";
  }
}


/***** output solution in given names of vertices and items *****/
#define TABLE_MAX 3500
#define TABNAME_MAX 32
void outputByTable(Param P, Tool T, Graph G, SeqForest &SeqF){  
  char *str, *foo, *bar;
  FILE *fp_out,*fp_v,*fp_i;
  OrderedSeqTree RB;  
  char *VT[TABLE_MAX], *IT[TABLE_MAX];
  VertexIDSeq seq;
  int x,p,all_n=0,all_items=0;
  
  if(P->vtable == NULL || P->itable==NULL)
    return;
  str = new char[LINE_MAX];
  foo = new char[LINE_MAX];
  bar = new char[LINE_MAX];
  
  // read vtable
  fp_v = open_file(P->vtable, "r");
  while(fgets(str, LINE_MAX-1, fp_v)!=NULL){
    if(str[0]=='#')
      continue;
    sscanf(str, "%d%s%s%s", &p, foo, bar, foo);
    VT[p] = new char[TABNAME_MAX];
    strcpy(VT[p], bar);
    all_n++;
    if(all_n==TABLE_MAX){
      fprintf(stderr, "%s contains too many vertices.\n", P->vtable);
      exit(EXIT_FAILURE);
    }
  }
  fclose(fp_v);
  
  // read itable
  fp_i = open_file(P->itable, "r");
  while(fgets(str, LINE_MAX-1, fp_i)!=NULL){
    if(str[0]=='#')
      continue;
    sscanf(str, "%d%s%s", &p, foo, bar);
    IT[p] = new char[TABNAME_MAX];
    strcpy(IT[p], bar);
    all_items++;
    if(all_items==TABLE_MAX){
      fprintf(stderr, "%s contains too many vertices.\n", P->itable);
      exit(EXIT_FAILURE);
    }
  }
  fclose(fp_i);

  // write the result
  fp_out = open_file(P->outname, "w");
  for(auto treeitr=SeqF.begin(); treeitr!=SeqF.end(); ++treeitr){
    SeqTree tree = *treeitr;
    for(auto itr=tree.begin(); itr!=tree.end(); ++itr){
      seq = itr->first;
      sort(seq.begin(), seq.end());
      RB[seq] = true;
    }
  }
  for(auto itr=RB.begin(); itr!=RB.end(); ++itr){
    ItemIDSeq itmseq;
    seq = itr->first;
    itmseq = getItemIntersection(T, G, seq);
    fprintf(fp_out, "--- %lu vertices with %lu items ---\n\n", seq.size(), itmseq.size());
    int r = 0;
    fprintf(fp_out, "Vertices:\t{");
    for(auto itr2=seq.begin(); itr2!=seq.end(); ++itr2){
      x = *itr2;
      fprintf(fp_out, " %s", VT[T->VMapInv[G->V[x]->id]]);
      r++;
      if(r%5==0)
	fprintf(fp_out, "\n");
    }
    r = 0;
    fprintf(fp_out, "}\nItems:\t{");
    for(auto itr2=itmseq.begin(); itr2!=itmseq.end(); ++itr2){
      x = *itr2;
      fprintf(fp_out, " %s", IT[T->IMapInv[x]]);
      r++;
      if(r%8==0)
	fprintf(fp_out, "\n");
    }
    fprintf(fp_out, "}\n\n\n");
  }
  
  fclose(fp_out);
  delete[] str;
  delete[] foo;
  delete[] bar;
}


/***** output distribution of connectors *****/
void outputDist(FILE *fp, Tool T, Graph G, TrieForest &TrieF){
  for(auto itr=TrieF.begin(); itr!=TrieF.end(); itr++){
    Trie trie = *itr;
    TrieNode u = trie->first_leaf;
    while(u!=NULL){
      fprintf(fp, "%d\t%d\t%g\n", u->seq.size(), getNumItems(T->B, T->d, G, u->seq),
	      getDensity(u->seq.size(), u->deg/2));
      u = u->next;
    }
  }
}
/***** output G to stdout *****/
void outputGraph(Tool T, Graph G){
  int i,x;
  for(x=0;x<G->n;x++){
    if(G->V[x]==NULL)
      continue;
    printf("%d --->",T->VMapInv[G->V[x]->id]);
    AdjList list;
    list = G->V[x]->A;
    while(list!=NULL){
      printf(" %d",T->VMapInv[list->v->id]);
      list = list->next;
    }
    printf("\titemset: {");
    for(i=0;i<G->items;i++)
      if(hasItem(T->B, T->d, G->V[x]->I, i))
	printf(" %d",T->IMapInv[i]);
    printf("}\n");
  }
  printf("\n");
}

/***** WP3 algorithm *****/
//This uses fewer arithmetic operations than any other known  
//implementation on machines with fast multiplication.
//It uses 12 arithmetic operations, one of which is a multiply.
static unsigned WP3(BitString n)
{
    BitString m1 = (~(BitString)0) / 3;
    BitString m2 = (~(BitString)0) / 5;
    BitString m4 = (~(BitString)0) / 17;
    BitString h01 = (~(BitString)0) / 255;

    n -= (n >> 1) & m1;             //put count of each 2 bits into those 2 bits
    n = (n & m2) + ((n >> 2) & m2); //put count of each 4 bits into those 4 bits 
    n = (n + (n >> 4)) & m4;        //put count of each 8 bits into those 8 bits 

    return (n * h01) >> (BIT_LENGTH-8);
    // returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
}


double getDensity(int n, int m){
  return ((double)m/(double)n) * 2.0 / (double)(n-1);
}

/***** trie related functions *****/
void initTrie(Trie T, int h){
  T->root = new struct _TrieNode;
  T->root->leaf = false;
  T->root->key = 0;
  T->first_leaf = NULL;
  T->h = h;
  T->size = 0;
}

bool memberTrie(Trie T, VertexIDSeq &seq){
  bool flag = true;
  TrieNode u = T->root;
  TrieEdge e;
  while(!u->leaf){
    // node
    if(u->E.find(seq[u->key]) == u->E.end()){
      flag = false;
      break;
    }
    // edge
    e = u->E[seq[u->key]];
    for(int h=e->head;h<=e->tail;h++)
      if(seq[h]!=e->slice[h-e->head]){
	flag = false;
	break;
      }
    if(!flag)
      break;
    u = e->child;
  }
  return flag;
}

bool insertTrie(Trie T, VertexIDSeq &seq, int value, int deg){
  bool exist = true;
  TrieNode u = T->root;  
  TrieEdge e;
  while(!u->leaf){
    // node    
    if(u->E.find(seq[u->key]) == u->E.end()){
      exist = false;
      u->E[seq[u->key]] = new struct _TrieEdge;
      e = u->E[seq[u->key]];
      e->parent = u;
      e->child = new struct _TrieNode;
      e->head = u->key;
      e->tail = T->h-1;
      for(int h=e->head;h<=e->tail;h++)
	e->slice.push_back(seq[h]);
      u = e->child;
      break;
    }
    // edge
    e = u->E[seq[u->key]];
    for(int h=e->head;h<=e->tail;h++){
      if(seq[h]==e->slice[h-e->head])
	continue;
      exist = false;
      TrieNode v = new struct _TrieNode;
      TrieNode w = new struct _TrieNode;
      TrieEdge _e = new struct _TrieEdge;
      TrieEdge __e = new struct _TrieEdge;
      // _e
      _e->parent = v;
      _e->child = e->child;
      _e->head = h;
      _e->tail = e->tail;
      for(int l=h;l<=e->tail;l++)
	_e->slice.push_back(e->slice[l-e->head]);
      // __e
      __e->parent = v;
      __e->child = w;
      __e->head = h;
      __e->tail = T->h-1;
      for(int l=h;l<T->h;l++)
	__e->slice.push_back(seq[l]);
      // v
      v->leaf = false;
      v->key = h;
      v->E[_e->slice[0]] = _e;
      v->E[seq[h]] = __e;      
      // e: parent and head do not change
      e->child = v;
      e->tail = h-1;
      // u
      u = w;
      break;
    }
    if(!exist)
      break;
    u = e->child;
  }
  if(exist==false){
    for(auto itr=seq.begin(); itr!=seq.end(); ++itr){
      u->seq.push_back(*itr);
    }
    u->next = T->first_leaf;
    u->leaf = true;
    u->deg = deg;
    T->first_leaf = u;
    T->size++;
  }
  u->value = value;
  return exist;
}

/***** output TrieForest *****/
void outputTrieForest(Param P, Tool T, Graph G, TrieForest &TrieF){
  for(auto itr=TrieF.rbegin(); itr!=TrieF.rend(); ++itr){
    SeqTree SeqT;
    if((*itr)->h < P->sigma)
      continue;
    TrieNode u = (*itr)->first_leaf;
    while(u!=NULL){
      SeqT[u->seq] = 1; 
      u = u->next;
    }
    outputSeqTree(T, G, SeqT);
  }
  cout << "\n";
}


#if 0
unsigned long long getNumNodes(TrieNode u){
  unsigned long long num=0;
  for(auto itr=u->c.begin(); itr!=u->c.end(); itr++)
    num += getNumNodes(itr->second);
  return 1+num;
}

void initTrie(Trie T, int h){
  T->root = new struct _TrieNode;
  T->first_leaf = NULL;
  T->h = h;
  T->size = 0;
}

bool memberTrie(Trie T, VertexIDSeq &seq){
  bool flag = true;
  TrieNode u;
  u = T->root;
  for(int h=0;h<T->h;h++){
    if(u->c.find(seq[h]) == u->c.end()){
      flag = false;
      break;
    }
    u = u->c[seq[h]];
  }
  return flag;
}

bool insertTrie(Trie T, VertexIDSeq &seq, int value){
  bool exist = true;
  TrieNode u;
  u = T->root;
  for(int h=0;h<T->h;h++){
    if(u->c.find(seq[h]) == u->c.end()){
      u->c[seq[h]] = new struct _TrieNode;
      exist = false;
    }
    u = u->c[seq[h]];
  }
  if(exist==false){
    for(auto itr=seq.begin(); itr!=seq.end(); ++itr){
      u->seq.push_back(*itr);
    }
    u->next = T->first_leaf;
    T->first_leaf = u;
    T->size++;
  }
  u->value = value;
  return exist;
}
#endif
