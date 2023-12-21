/********** define.h **********/

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

#define BIT_LENGTH (sizeof(BitString)*8)
typedef unsigned long long BitString;
typedef BitString *ItemStr; 

typedef unsigned long long HugePositive;

typedef vector<int> VertexIDSeq;
typedef vector<int> ItemIDSeq;
typedef vector<int> IntvIDSeq;

// other typedefs and structures
typedef struct _AdjList *AdjList;
typedef struct _Vertex *Vertex;
typedef struct _Component *Component;
typedef struct _Edge *Edge;
typedef struct _Graph *Graph;
typedef struct _Tool *Tool;
typedef struct _BFSTool *BFSTool;
typedef struct _Param *Param;

struct _Vertex{
  int id;             // vertex id
  AdjList A;          // adjacency list; pointer to list element
  int deg;            // degree of vertex
  ItemStr I;          // bit string that represents item set
  IntvIDSeq *Itv_ptr; // pointers to IDs of bit strings that contain at least one item
  int items;          // number of one's in I
};

struct _Component{
  VertexIDSeq seq; // vertex ids in the component
  int minitem;     // min item of the component
  ItemStr I;       // bit string that represents the common item set
  IntvIDSeq Itv;   // IDs of bit strings that contain at least one item
};

struct _AdjList{
  Vertex v;      // pointer to vertex
  double w;      // edge weight
  AdjList next;  // next list element
};

struct _Graph{
  Vertex *V;     // array of vertices
  int n;         // number of vertices
  int m;         // number of edges
  int items;     // number of items
  double item_density; // avg of v->items / items
};

struct _Tool{
  BitString B[BIT_LENGTH]; // base vector (B[j][j]=1 and B[j][k]=0 for k\ne j)
  BitString H[BIT_LENGTH]; // (head)1000...0(tail), 1100...0, ..., 1111...1
  BitString T[BIT_LENGTH]; // (head)1111...1(tail), 0111...1, ..., 000...01
  int d;                   // = q/BIT_LENGTH + 1  
  unordered_map<int, int> VMap;    // vertex name to vertex id
  unordered_map<int, int> VMapInv; // vertex id to vertex name
  unordered_map<int, int> IMap;    // item name to item id
  unordered_map<int, int> IMapInv; // item id to item name
  int reduced_vertices;
  int reduced_edges;
  VertexIDSeq seq_all;
  vector<VertexIDSeq> seq_item;
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
#define INI_tlim -1.0
#define INI_reduce true
#define INI_outname "out.txt"
#define INI_distname "dist.txt"
#define INI_ramub -1

struct _Param{
  char *ptn_file;
  char *grh_file;
  int theta;
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
