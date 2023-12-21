/********** define.h **********/


// macros
#ifndef LINE_MAX
#define LINE_MAX 4096
#endif 

#define LEN_FILENAME 256

#define MARKER_MAX INT_MAX-1

#define EVERY_ITEM -1
#define BFS_WHITE 0
#define BFS_GRAY 1
#define BFS_BLACK 2
#define Undef -1

#define Smaller -1
#define Equiv 0
#define Larger 1

#define BIT_LENGTH sizeof(BitString)*8
typedef unsigned long long BitString;
typedef BitString *ItemStr; 

/*
  // for future update (2018/6/12)
  typedef int16_t VertexID;
  typedef int16_t ItemID;
*/

typedef vector<int> VertexIDSeq;
typedef vector<int> ItemIDSeq;

struct ordered_struct{
  bool operator() (const VertexIDSeq &left, const VertexIDSeq &right) const
  {
    int i;
    if(left.size() > right.size())
      return false;
    if(left.size() < right.size())
      return true;
    for(i=0;i<left.size();i++)
      if(left[i]<right[i])
        return true;
      else if(left[i]>right[i])
        return false;
    return false;
  }
};

struct unordered_struct{
  bool operator() (const VertexIDSeq &left, const VertexIDSeq &right) const
  {
    int i,size;
    size = std::min(left.size(), right.size());
    for(i=0;i<size;i++)
      if(left[i]<right[i])
        return true;
      else if(left[i]>right[i])
        return false;
    return false;
  }
};

typedef map<VertexIDSeq, int, ordered_struct> OrderedSeqTree;
typedef map<VertexIDSeq, int, unordered_struct> SeqTree;
typedef vector<SeqTree> SeqForest;

/********** Trie related data structures **********/
typedef struct _TrieNode *TrieNode;
typedef struct _TrieEdge *TrieEdge;
typedef struct _Trie *Trie; 
typedef unordered_map<int, TrieEdge> ChildEdge;
typedef vector<Trie> TrieForest;

struct _TrieNode{
  bool leaf;
  
  /*** used only when leaf=true ***/
  VertexIDSeq seq;
  TrieNode next;
  int value;
  int deg;
  
  /*** used only when leaf=false ***/
  int key;
  ChildEdge E;
};

struct _TrieEdge{
  TrieNode parent;
  TrieNode child;
  int head;
  int tail;
  VertexIDSeq slice;
};

struct _Trie{
  TrieNode root;
  TrieNode first_leaf;
  int h;
  unsigned long long size;
};
/**************************/

// other typedefs and structures
typedef struct _AdjList *AdjList;
typedef struct _Vertex *Vertex;
typedef struct _Edge *Edge;
typedef struct _Graph *Graph;
typedef struct _Tool *Tool;
typedef struct _Param *Param;

struct _Vertex{
  int id;        // vertex id
  ItemStr I;     // bit string that represents item set
  AdjList A;     // adjacency list; pointer to list element
  int mark;      // mark for BFS
  int bfs_color; // color for BFS
  int bfs_conn;  // id of connected component in BFS
  int items;     // number of one's in I
};

struct _Edge{
  Vertex V[2];   // pointers to extreme points
  double w;      // edge weight
};

struct _AdjList{
  Vertex v;      // pointer to vertex
  Edge e;        // pointer to edge
  AdjList next;  // next list element
};

struct _Graph{
  Vertex *V;     // array of vertices
  Edge *E;       // array of edges
  int n;         // number of vertices
  int m;         // number of edges
  int items;     // number of items
  double item_density; // avg of v->items / items
};

struct _Tool{
  BitString B[BIT_LENGTH];
  int d;                   // = q/BIT_LENGTH + 1
  unsigned int marker;     // marker for BFS in {1,...,marker_max-1}
  unsigned int marker_max; // will be initialized to MARKER_MAX
  unordered_map<int, int> VMap;
  unordered_map<int, int> VMapInv;
  unordered_map<int, int> IMap;
  unordered_map<int, int> IMapInv;
  int reduced_vertices;
  int reduced_edges;
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
