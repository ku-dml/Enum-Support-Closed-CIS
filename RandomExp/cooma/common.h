/***** common.h *****/
void checkArgs(int argc, char *argv[]);
void readArgs(Param P, Graph G, int argc, char *argv[]);
void preprocess(Param P, Tool T, Graph G);
void readPattern(Param P, Tool T, Graph G);
void readGraph(Param P, Tool T, Graph G);
int countNumEdges(Graph G);
void addItem(BitString B[BIT_LENGTH], int d, ItemStr I, int i);
bool hasItem(BitString B[BIT_LENGTH], int d, ItemStr I, int i);
bool isUnionSmall(BitString B[BIT_LENGTH], int d, Graph G, VertexIDSeq seq, int theta);
bool isIntersectionSmall(BitString B[BIT_LENGTH], int d, Graph G, VertexIDSeq seq, int theta);
int getNumItems(BitString B[BIT_LENGTH], int d, Graph G, VertexIDSeq seq);

ItemIDSeq getItemIntersection(Tool T, Graph G, VertexIDSeq seq);
VertexIDSeq getSubseqByItem(Tool T, Graph G, int i, VertexIDSeq seq);
VertexIDSeq getIntersection(Tool T, Graph G, VertexIDSeq seq);
Graph getSubgraph(Graph G, VertexIDSeq seq);
void deleteSubgraph(Graph G);
int markSet(Tool T, Graph U, Graph G, VertexIDSeq seq);
bool increaseMarker(Tool T, Graph U, Graph G);
int markedBFS(Tool T, Graph G, int marker, Vertex v, int conn);
SeqTree computeSeqTree(Tool T, Graph G, VertexIDSeq seq, int marker, int i);
void outputSet(Tool T, Graph G, VertexIDSeq seq);
void outputSeqTree(Tool T, Graph G, SeqTree &SeqT);
void outputSeqForest(Tool T, Graph G, SeqForest &SeqF);
void outputInStdform(Tool T, Graph G, SeqForest &SeqF);
void outputByTable(Param P, Tool T, Graph G, SeqForest &SeqF);
void outputDist(FILE *fp, Tool T, Graph G, TrieForest &TrieF);
void outputGraph(Tool T, Graph G);
static unsigned WP3(BitString n);
double getDensity(int n, int m);

/***** trie related functions *****/
void initTrie(Trie T, int h);
bool memberTrie(Trie T, VertexIDSeq &seq);
bool insertTrie(Trie T, VertexIDSeq &seq, int value, int deg);
void outputTrieForest(Param P, Tool T, Graph G, TrieForest &TrieF);
