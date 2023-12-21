/***** common.h *****/
void checkArgs(int argc, char *argv[]);
void readArgs(Param P, Graph G, int argc, char *argv[]);
void preprocess(Param P, Tool T, Graph G);
void readPattern(Param P, Tool T, Graph G);
void readGraph(Param P, Tool T, Graph G);
void initTool(Tool T, Graph G);
int countNumEdges(Graph G);

void addItem(BitString B[BIT_LENGTH], int d, ItemStr I, int i);
bool hasItem(BitString B[BIT_LENGTH], ItemStr I, int i);
VertexIDSeq getSubseqByItem(Tool T, Graph G, int i, VertexIDSeq &seq);

int getCmaxByItem(Tool T, Graph G, BFSTool B, Component Parent, int j, vector<Component> &Cmax);
int prepareBFSByItem(Tool T, BFSTool B, Graph G, unsigned int *marker_ptr, VertexIDSeq &seq, int j);


int getCmax(Tool T, Graph G, BFSTool B, VertexIDSeq &seq,
	    vector<Component> &Cmax);
unsigned int prepareBFS(BFSTool B, int n, VertexIDSeq &seq);
void markedBFS(BFSTool B, unsigned int marker, Vertex v, VertexIDSeq &S);
bool isLarger(BFSTool B, unsigned int marker, Vertex v, int num);
bool isMinItemK(Tool T, Graph G, Component C, int k);
bool isCommonItemEmpty(Tool T, Graph G, Component C);
int getFirstBit(BitString H[BIT_LENGTH], BitString T[BIT_LENGTH], BitString r, int left, int right);
void getItemIntersection(Tool T, Graph G, Component C);

int getDiffMinItem(Tool T, Graph G, int k, Component Parent, Component S);
int getDiffMinItem(Tool T, Graph G, int k, ItemStr T_I, ItemStr S_I);

bool isParent(Tool T, Graph G, BFSTool B, Component S, Component Parent); 
bool isParent(Tool T, Graph G, BFSTool B, Component S, Component Parent, vector<VertexIDSeq>&Pseq); 

void getSubgraph(Tool T, Graph G, Graph S, int k);
void delSubgraph(Graph S);

void outputGraph(Tool T, Graph G);
static unsigned WP3(BitString n);
double getDensity(int n, int m);

