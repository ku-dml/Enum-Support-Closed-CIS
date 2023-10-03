/********** define.h **********/
#ifndef DEFINE_H
#define DEFINE_H


#include <bitset>
#include <set>
#include <unordered_map>
#include <queue>

// macros
#ifndef LINE_MAX
#define LINE_MAX 65536
#endif 

#define LEN_FILENAME 256

#define MARKER_MAX INT_MAX-1

#define BFS_WHITE 0
#define BFS_GRAY 1
#define BFS_BLACK 2
#define Undef -1 // this must be -1 

#define Smaller -1
#define Equiv 0
#define Larger 1
#define ITEM_SIZE 500 //101
#define VERTEX_SIZE 100001

#define BIT_LENGTH (sizeof(BitString)*8)

using namespace std;

typedef unsigned long long BitString;
typedef BitString *ItemStr;



typedef unsigned long long HugePositive;

typedef vector<int> VertexIDSeq;
typedef bitset<VERTEX_SIZE> vertexset;
typedef vector<int> ItemIDSeq;
typedef bitset<ITEM_SIZE> itemset;
typedef vector<int> IntvIDSeq;



// other typedefs and structures
typedef struct _AdjList *AdjList;
typedef struct _Vertex *Vertex;
typedef struct _Edge *Edge;
typedef struct _Graph *Graph;
typedef struct _Tool *Tool;
typedef struct _BFSTool *BFSTool;
typedef struct _Param *Param;

typedef struct _OwnStack *OwnStack;

typedef _OwnStack BanList;
typedef _OwnStack Solution;

struct _Vertex{
  int id;             // vertex id
  VertexIDSeq A;          // adjacency list;List of vertex ID
  int deg;            // degree of vertex
  itemset *I;          // bit string that represents item set
  IntvIDSeq *Itv_ptr; // pointers to IDs of bit strings that contain at least one item
  int items;          // number of one's in I
};


struct _Graph{
  Vertex *V;     // array of vertices
  int n;         // number of vertices
  int m;         // number of edges
  int items;     // number of items
  double item_density; // avg of v->items / items
};

struct _Tool{
  unordered_map<int, int> VMap;    // vertex name to vertex id
  unordered_map<int, int> VMapInv; // vertex id to vertex name
  unordered_map<int, int> IMap;    // item name to item id
  unordered_map<int, int> IMapInv; // item id to item name
  unordered_map<int, int> PMap;    // population name to population id
  unordered_map<int, int> PMapInv; // population id to population name
  int reduced_vertices;
  int reduced_edges;
  VertexIDSeq seq_all;
  vector<VertexIDSeq> seq_item;
  itemset phenotype; 
  vector<itemset> population;
  int numPop ;
  vector<int> n;	//nj in the paper
  vector<int> n1;
  vector<int> n2;
};

struct _BFSTool{
  unsigned int marker;     // marker for BFS in {1,...,marker_max-1}
  unsigned int marker_max; // will be initialized to MARKER_MAX
  vector<unsigned int> mark;        // mark for BFS
  vector<unsigned int> bfs_color;   // color for BFS
};

#define INI_seed 1
#define INI_sigma 1
#define INI_delta 0.0
#define INI_turnWidth 1;
#define INI_tlim -1.0
#define INI_reduce true
#define INI_outname "out.txt"
#define INI_common_outname "results.csv"
#define INI_distname "dist.txt"
#define INI_ramub -1


struct _Param{
  char *ptn_file;
  char *grh_file;
  char *phenotype_file;
  char *population_file; 
  int theta;
  double alpha;
  int turnWidth;
  int seed;
  int sigma;
  double delta;
  double tlim; // time limit
  bool reduce; // whether reduction is invoked
  char *vtable;
  char *itable;
  char *outname;
  char *common_outname;
  char *distname;
  int ramub;   // upper bound on RAM
};

struct _OwnStack{
  VertexIDSeq seq;
  vertexset bit;
  
  void pop_back(){
    if (seq.back() != -1) {
      bit.reset(seq.back());
    }
    seq.pop_back();	
  }
  
  void push_back(int key){
    seq.push_back(key);
    if (key != -1) {
      bit.set(key);
    }
  }

};
#endif
