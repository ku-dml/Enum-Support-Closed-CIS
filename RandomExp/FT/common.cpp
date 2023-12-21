/********** common.cpp **********/

/***** includes *****/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

// C++ libraries
#include <iostream>
#include <unordered_map>
#include <queue>
#include <vector>
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
  P->tlim = INI_tlim; // negative means infinity
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
      fprintf(stderr, "error: <%s> is illegal parameter.\n", arg);
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
  int b,i,j,items;

  str = new char[LINE_MAX]; 
  str_id = new char[LINE_MAX];
  str_item = new char[LINE_MAX];

  /* first scan of pattern file */
  G->n = 0;
  //G->m = 0;
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
  
  
  /* initialize T */
  T->B[0] = 1;
  for(b=1;b<BIT_LENGTH;b++)
    T->B[b] = T->B[b-1] << 1;
  T->H[0] = 1;
  for(b=1;b<BIT_LENGTH;b++)
    T->H[b] = 1 + (T->H[b-1] << 1);
  T->T[0] = (~(BitString)0);
  for(b=1;b<BIT_LENGTH;b++)
    T->T[b] = T->T[b-1] << 1;  
  T->d = G->items/BIT_LENGTH + 1;  
  T->reduced_vertices = 0;
  T->reduced_edges = 0;
  for(j=0;j<G->items;j++){
    VertexIDSeq empty;
    T->seq_item.push_back(empty);
  }
  
  delete[] str;
  delete[] str_id;
  delete[] str_item;
}

void initTool(Tool T, Graph G){
  int i,j;
  for(i=0;i<G->n;i++)
    if(G->V[i]!=NULL)
      T->seq_all.push_back(i);
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
    G->V[x]->A = NULL;
    G->V[x]->deg = 0;
    G->V[x]->Itv_ptr = new IntvIDSeq[1];
    G->V[x]->I = new BitString[T->d];
    for(p=0;p<T->d;p++)
      G->V[x]->I[p] = 0;
    if(strcmp(str_item, "")==Equiv || strcmp(str_item, "!")==Equiv)
      G->V[x]->items = 0;
    else{
      s = split(str_item, ',');
      G->V[x]->items = getNumChar(str_item, ',')+1;
      item_density += (double)(G->V[x]->items);
      for(j=0;j<G->V[x]->items;j++){
	i = T->IMap[atoi(s[j])];
	addItem(T->B, T->d, G->V[x]->I, i);
	T->seq_item[i].push_back(x);
	free(s[j]);
      }
      free(s);
      for(p=0;p<T->d;p++)
	if(G->V[x]->I[p])
	  G->V[x]->Itv_ptr->push_back(p);
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
  int a,b,x,y,line=0;
  
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

    AdjList list;
    bool added = true;
    
    list = new struct _AdjList;
    list->v = G->V[y];
    list->w = w;

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
	  delete list;
	  added = false;
	}
	else{
	  list->next = prev->next;
	  prev->next = list;
	}
      }
    }
    if(added)
      (G->V[x]->deg)++;
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
bool hasItem(BitString B[BIT_LENGTH], ItemStr I, int i){
  int p;
  p = i/BIT_LENGTH;
  if(B[i%BIT_LENGTH]&I[p])
    return true;
  return false;
}


/***** get subsequence of seq such that 
       the vertices have item i *****/
VertexIDSeq getSubseqByItem(Tool T, Graph G, int i, VertexIDSeq &seq){
  VertexIDSeq subseq;
  int x,k,size;
  size=seq.size();
  for(k=0;k<size;k++){
    x = seq[k];
    if(G->V[x]!=NULL && hasItem(T->B, G->V[x]->I, i))
      subseq.push_back(x);
  }
  return subseq;
}

/***** compute collection of maximal components: faster version *****/
int getCmaxByItem(Tool T, Graph G, BFSTool B, Component Parent, int j, vector<Component> &Cmax){
  unsigned int marker;
  int size,numComp=0;
  prepareBFSByItem(T, B, G, &marker, Parent->seq, j);
  size = Parent->seq.size();
  for(int i=0;i<size;i++)
    if(B->bfs_color[Parent->seq[i]] == BFS_WHITE &&
       B->mark[Parent->seq[i]] == marker){
      Component C;
      C = new struct _Component;
      markedBFS(B, marker, G->V[Parent->seq[i]], C->seq);
      C->I = NULL;
      Cmax.push_back(C);
      numComp++;
    }
  return numComp;
}

int prepareBFSByItem(Tool T, BFSTool B, Graph G, unsigned int *marker_ptr, VertexIDSeq &seq, int j){
  unsigned int marker;
  int size,num=0;  
  *marker_ptr = B->marker;
  B->marker++;
  // reset the marks of all vertices if marker=marker_max
  if(B->marker>B->marker_max){
    for(int x=0;x<G->n;x++)
      B->mark[x] = 0;
    B->marker = 1;
  }
  size=seq.size();
  for(int i=0;i<size;i++)
    if(hasItem(T->B, G->V[seq[i]]->I, j)){
      B->mark[seq[i]] = *marker_ptr;
      B->bfs_color[seq[i]] = BFS_WHITE;
      num++;
    }
  return num;
}


/***** compute collection of maximal components *****/
int getCmax(Tool T, Graph G, BFSTool B, VertexIDSeq &seq, vector<Component> &Cmax){
  unsigned int marker;
  int size,numComp=0;
  marker = prepareBFS(B, G->n, seq);
  size = seq.size();
  for(int i=0;i<size;i++)
    if(B->bfs_color[seq[i]] == BFS_WHITE){
      Component C;
      C = new struct _Component;
      markedBFS(B, marker, G->V[seq[i]], C->seq);
      C->I = NULL;
      Cmax.push_back(C);
      numComp++;
    }
  return numComp;
}

/***** prepare BFS *****/
unsigned int prepareBFS(BFSTool B, int n, VertexIDSeq &seq){
  unsigned int marker;
  int size;  
  marker = B->marker;
  B->marker++;
  // reset the marks of all vertices if marker=marker_max
  if(B->marker>B->marker_max){
    for(int x=0;x<n;x++)
      B->mark[x] = 0;
    B->marker = 1;
  }
  size=seq.size();
  for(int i=0;i<size;i++){
    B->mark[seq[i]] = marker;
    B->bfs_color[seq[i]] = BFS_WHITE;
  }
  return marker;
}


/***** BFS that starts from vertex v in S
       it returns minitem in C *****/
void markedBFS(BFSTool B, unsigned int marker, Vertex v, VertexIDSeq &S){
  AdjList A;
  Vertex u;
  queue<Vertex> Q;
  Q.push(v);
  B->bfs_color[v->id] = BFS_GRAY;
  while(Q.empty()==false){
    u = Q.front();
    B->bfs_color[u->id] = BFS_BLACK;
    Q.pop();
    S.push_back(u->id);
    A = u->A;
    while(A!=NULL){
      if(B->mark[A->v->id]==marker && B->bfs_color[A->v->id]==BFS_WHITE){
	Q.push(A->v);
	B->bfs_color[A->v->id] = BFS_GRAY;
      }
      A = A->next;
    }
  }
}

/***** identify whether v belongs to a component whose size > num *****/
bool isLarger(BFSTool B, unsigned int marker, Vertex v, int num){
  AdjList A;
  Vertex u;
  queue<Vertex> Q;
  int num_black = 0;
  Q.push(v);
  B->bfs_color[v->id] = BFS_GRAY;
  while(Q.empty()==false){
    u = Q.front();
    B->bfs_color[u->id] = BFS_BLACK;
    Q.pop();
    num_black++;
    if(num_black>num)
      return true;
    A = u->A;
    while(A!=NULL){
      if(B->mark[A->v->id]==marker && B->bfs_color[A->v->id]==BFS_WHITE){
	Q.push(A->v);
	B->bfs_color[A->v->id] = BFS_GRAY;
      }
      A = A->next;
    }
  }
  return false;
}


/***** decide whether min(C->seq's itemset)=k;
       when true, C->Itv is determined correctly.
       It is assumed that all vertices in C have item k *****/
bool isMinItemK(Tool T, Graph G, Component C, int k){
  BitString b;
  int p_of_k,size;
  IntvIDSeq::iterator itr_itv;
  
  p_of_k = k/BIT_LENGTH;
  size = C->seq.size();
 
  // determine the common item strings for intervals < p_of_k
  for(itr_itv=G->V[C->seq[0]]->Itv_ptr->begin();
      itr_itv!=G->V[C->seq[0]]->Itv_ptr->end();
      itr_itv++){
    int p=*itr_itv;
    if(p>=p_of_k)
      break;
    b = ~(BitString)0;
    for(int i=0;i<size;i++){
      b = b & G->V[C->seq[i]]->I[p];
      if(b==0)
	break;
    }
    if(b!=0)
      return false;
  }
  
  // decide whether minitem is k or not
  b = ~(BitString)0;
  for(int i=0;i<size;i++)
    b = b & G->V[C->seq[i]]->I[p_of_k];
  int h;
  h = getFirstBit(T->H, T->T, b, 0, BIT_LENGTH-1);
  if(p_of_k*BIT_LENGTH+h!=k)
    return false;

  // determine the common item strings for intervals > p_of_k
#if 1
  C->Itv.push_back(p_of_k);
  for(;itr_itv!=G->V[C->seq[0]]->Itv_ptr->end();itr_itv++){
    int p=*itr_itv;
    if(p==p_of_k)
      continue;
    b = ~(BitString)0;
    for(int i=0;i<size;i++){
      b = b & G->V[C->seq[i]]->I[p];
      if(b==0)
	break;
    }
    if(b!=0)
      C->Itv.push_back(p);
  }
#else
  int *P;
  P = new int[T->d];
  for(int p=p_of_k;p<T->d;p++)
    P[p] = 0;

    
  for(int i=0;i<size;i++){
    int x = C->seq[i];
    for(auto itr_itv=G->V[x]->Itv_ptr->rbegin();
	itr_itv!=G->V[x]->Itv_ptr->rend();
	itr_itv++){
      int p = *itr_itv;
      if(p<p_of_k)
	break;
      P[p]++;
    }
  }

  // construct C->Itv; this shouldn't contain intervals 0 to p_of_k-1
  for(int p=p_of_k;p<T->d;p++)
    if(P[p]==size)
      C->Itv.push_back(p);
  delete[] P;
#endif
  
  C->minitem = k;
  return true;
}
  
bool isCommonItemEmpty(Tool T, Graph G, Component C){
  int size = C->seq.size();
  for(int p=0;p<T->d;p++){
    BitString b = ~(BitString)0;
    for(int i=0;i<size;i++){
      b = b & G->V[C->seq[i]]->I[p];
      if(b==0)
	break;
    }
    if(b!=0)
      return false;
  }
  return true;
}


int getFirstBit(BitString H[BIT_LENGTH], BitString T[BIT_LENGTH], BitString r, int left, int right){
  int mid;
  r = r&T[left];
  if(r==0)
    return Undef;
  while(1){
    mid = (left+right)/2;
    if(r&H[mid]){
      if(mid==left || !(r&H[mid-1]))
	break;
      right = mid-1;
    }
    else
      left = mid+1;
  }
  return mid;
}


void getItemIntersection(Tool T, Graph G, Component C){
  IntvIDSeq::iterator itr = C->Itv.begin();
  C->I = new BitString[T->d];
  for(int p=0;p<T->d;p++)
    if(itr==C->Itv.end())
      C->I[p] = 0;
    else{
      if(p==*itr){
	C->I[p] = (~(BitString)0);
	for(auto v_itr=C->seq.begin();v_itr!=C->seq.end();v_itr++)
	  C->I[p] = C->I[p] & G->V[*v_itr]->I[p];
	int l=0,j;
	/*
	while(1){
	  l = getFirstBit(T->H, T->T, C->I[p],
			  l, BIT_LENGTH-1);
	  if(l==Undef)
	    break;
	  j = p*BIT_LENGTH+l;
	  l++;
	  if(l==BIT_LENGTH)
	    break;
	}
	*/
	itr++;
      }
      else if(p<*itr)
	C->I[p] = 0;
      else{
	fprintf(stderr, "error: something wrong happens at getItemIntersection.\n");
	exit(1);
      }
    }
}

int getDiffMinItem(Tool T, Graph G, int k, Component Parent, Component S){
  BitString b,x;
  int p,p0,q,i,a=Undef;
  p0 = (k+1)/BIT_LENGTH;
  q = (k+1)%BIT_LENGTH;
  for(auto itr=S->Itv.begin();itr!=S->Itv.end();itr++){
    p = *itr;
    x = S->I[p] & Parent->I[p];
    b = S->I[p] ^ x;// by this, b = S->I[p]-Parent->I[p]
    if(p==p0)
      b = b & T->T[q];
    if(b){      
      i = getFirstBit(T->H, T->T, b, 0, BIT_LENGTH-1);
      a = p*BIT_LENGTH + i;
      break;
    }
  }
  return a;

}

int getDiffMinItem(Tool T, Graph G, int k, ItemStr T_I, ItemStr S_I){
  BitString b,x;
  int p0,q,i,a=Undef;
  p0 = (k+1)/BIT_LENGTH;
  q = (k+1)%BIT_LENGTH;
  for(int p=p0;p<T->d;p++){
    x = S_I[p] & T_I[p];
    b = S_I[p] ^ x;       // by this, b = S_I[p]-T_I[p]
    if(p==p0)
      b = b & T->T[q];
    if(b){      
      i = getFirstBit(T->H, T->T, b, 0, BIT_LENGTH-1);
      a = p*BIT_LENGTH + i;
      break;
    }
  }
  return a;
}

#if 1
bool isParent(Tool T, Graph G, BFSTool B, Component S, Component Parent, vector<VertexIDSeq>&Pseq){
  bool flag=true;
  BitString b;
  int j,l,numS,t=0;
  VertexIDSeq seq,V;
  unsigned int marker;

  numS = S->seq.size();
 
  for(auto itr=S->Itv.begin(); itr!=S->Itv.end(); itr++){
    int p = *itr;
    b = S->I[p];
    l = 0;
    while(1){
      l = getFirstBit(T->H, T->T, b, l, BIT_LENGTH-1);
      if(l==Undef)
        break;
      j = p*BIT_LENGTH+l;
      if(j==S->minitem){
	l++;
	if(l==BIT_LENGTH)
	  break;
	continue;
      }
      if(hasItem(T->B, Parent->I, j))
	t++;
      else{
	int numT;
	numT = prepareBFSByItem(T, B, G, &marker, Pseq[t], j);
	if(numT!=numS)
	  flag = !isLarger(B, marker, G->V[S->seq[0]], numS);
      }
      
      
      if(flag==false)
	break;
      l++;
      if(l==BIT_LENGTH)
        break;
    }
    if(flag==false)
      break;
  }

  return flag;
}
#else
bool isParent(Tool T, Graph G, BFSTool B, Component S, Component Parent){
  bool flag=true,first=true,has;
  BitString b;
  int j,l,numS;
  VertexIDSeq seq,V;
  unsigned int marker;

  numS = S->seq.size();
 
  for(auto itr=S->Itv.begin(); itr!=S->Itv.end(); itr++){
    int p = *itr;
    b = S->I[p];
    l = 0;
    while(1){
      l = getFirstBit(T->H, T->T, b, l, BIT_LENGTH-1);
      if(l==Undef)
        break;
      j = p*BIT_LENGTH+l;
      if(j==S->minitem){
	l++;
	if(l==BIT_LENGTH)
	  break;
	continue;
      }
      has = hasItem(T->B, Parent->I, j);
      if(has){
	if(first){
	  seq = getSubseqByItem(T, G, j, T->seq_item[S->minitem]);
	  first = false;
	}
	else
	  seq = getSubseqByItem(T, G, j, seq);
      }
      else{
	int numT;
	if(first)
	  numT = prepareBFSByItem(T, B, G, &marker, T->seq_item[S->minitem], j);
	else
	  numT = prepareBFSByItem(T, B, G, &marker, seq, j);
	if(numT!=numS)
	  flag = !isLarger(B, marker, G->V[S->seq[0]], numS);
      }
      
      
      if(flag==false)
	break;
      l++;
      if(l==BIT_LENGTH)
        break;
    }
    if(flag==false)
      break;
  }

  return flag;
}
#endif

/***** get the subgraph of G induced by item k *****/
void getSubgraph(Tool T, Graph G, Graph S, int k){
  S->V = new Vertex[G->n];
  for(int i=0;i<G->n;i++){
    if(G->V[i]==NULL){
      S->V[i] = NULL;
      continue;
    }
    if(!(hasItem(T->B, G->V[i]->I, k))){
      S->V[i] = NULL;
      continue;
    }
    S->V[i] = new struct _Vertex;
    S->V[i]->id = i;
    S->V[i]->A = NULL; // will be constructed later
    S->V[i]->deg = 0;
    S->V[i]->I = G->V[i]->I;
    S->V[i]->Itv_ptr = G->V[i]->Itv_ptr;
    S->V[i]->items = G->V[i]->items;
  }
  for(int i=0;i<G->n;i++){
    if(S->V[i]==NULL)
      continue;
    AdjList list = G->V[i]->A;
    AdjList prev = NULL;
    while(list!=NULL){
      if(S->V[list->v->id]==NULL){
	list = list->next;
	continue;
      }
      AdjList new_list = new struct _AdjList;
      new_list->v = S->V[list->v->id];
      new_list->w = list->w;
      new_list->next = NULL;
      if(prev==NULL)
	S->V[i]->A = new_list;
      else
	prev->next = new_list;
      prev = new_list;
      list = list->next;
    }
  }
  S->n = G->n;
  S->m = G->m;
  S->items = G->items;
  S->item_density = G->item_density;
}


/***** delete inside of subgraph S *****/
void delSubgraph(Graph S){
  for(int i=0;i<S->n;i++){
    if(S->V[i]==NULL)
      continue;
    AdjList list,prev;
    list = S->V[i]->A;
    while(list!=NULL){
      prev = list;
      list = list->next;
      delete prev;
    }
    delete S->V[i];
  }
  delete[] S->V;
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
      if(hasItem(T->B, G->V[x]->I, i))
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

