/********** define.h **********/


#include <bitset>
#include <set>
#include <map>
#include <queue>
#include <unordered_map>

// macros
#ifndef LINE_MAX
#define LINE_MAX 65536
#endif

#ifndef DEFINE_H
#define DEFINE_H

#define LEN_FILENAME 256

#define MARKER_MAX INT_MAX-1

#define BFS_WHITE 0
#define BFS_GRAY 1
#define BFS_BLACK 2
#define Undef -1 // this must be -1 

#define Smaller -1
#define Equiv 0
#define Larger 1
#define ITEM_SIZE 1000
#define VERTEX_SIZE 20000

#define BIT_LENGTH (sizeof(BitString)*8)
typedef unsigned long long BitString;
typedef BitString *ItemStr;



typedef unsigned long long HugePositive;

typedef std::vector<int> VertexIDSeq;
typedef std::bitset<VERTEX_SIZE> vertexset;
typedef std::vector<int> ItemIDSeq;
typedef std::bitset<ITEM_SIZE> itemset;
typedef std::vector<int> IntvIDSeq;

// other typedefs and structures
typedef struct _AdjList *AdjList;
typedef struct _Vertex *Vertex;
typedef struct _Edge *Edge;
typedef struct _Graph *Graph;
typedef struct _Tool *Tool;
typedef struct _BFSTool *BFSTool;
typedef struct _Param *Param;

typedef struct _OwnStack *OwnStack;

struct _Vertex{
  int id;             // vertex id
  VertexIDSeq A;          // adjacency list;List of vertex ID
  int deg;            // degree of vertex
  std::bitset<ITEM_SIZE> *I;          // bit string that represents item set
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
  std::unordered_map<int, int> VMap;    // vertex name to vertex id
  std::unordered_map<int, int> VMapInv; // vertex id to vertex name
  std::unordered_map<int, int> IMap;    // item name to item id
  std::unordered_map<int, int> IMapInv; // item id to item name
  std::unordered_map<int, int> PMap;    // population name to population id
  std::unordered_map<int, int> PMapInv; // population id to population name
  int reduced_vertices;
  int reduced_edges;
  VertexIDSeq seq_all;
  std::vector<VertexIDSeq> seq_item;
  itemset phenotype; 
  std::vector<itemset> population;
  int numPop ;
  std::vector<int> n;	//nj in the paper
  std::vector<int> n1;
  std::vector<int> n2;

  OwnStack C;//Result
  OwnStack B;//check
};

struct _BFSTool{
  unsigned int marker;     // marker for BFS in {1,...,marker_max-1}
  unsigned int marker_max; // will be initialized to MARKER_MAX
  std::vector<unsigned int> mark;        // mark for BFS
  std::vector<unsigned int> bfs_color;   // color for BFS
};
  
#define INI_seed 1
#define INI_sigma 1
#define INI_delta 0.0
#define INI_tlim -1.0
#define INI_reduce true
#define INI_outname "out.txt"
#define INI_distname "dist.txt"
#define INI_ramub -1

struct _Param{
  char *ptn_file;
  char *grh_file;
  char *phenotype_file;
  char *population_file; 
  int theta;
  double alpha;
  int seed;
  int sigma;
  double delta;
  double tlim; // time limit
  bool reduce; // whether reduction is invoked
  char *vtable;
  char *itable;
  char *outname;
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