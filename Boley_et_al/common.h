/***** common.h *****/
#include "define.h"
#include <map>

void checkArgs(int argc, char *argv[]);
void readArgs(Param P, Graph G, int argc, char *argv[]);
void preprocess(Param P, Tool T, Graph G);
void readPattern(Param P, Tool T, Graph G);
void readPhenotype(Param P, Tool T);
void readPopulation(Param P, Tool T);
void readGraph(Param P, Tool T, Graph G);
void computeStatisticVal(Tool T);
void initTool(Tool T, Graph G);
int countNumEdges(Graph G);

unsigned int prepareBFS(BFSTool B, int n, VertexIDSeq &seq);

int getClosure(Param P, Tool T, Graph G, BFSTool B, Solution *S, BanList *Ban, Vertex V, itemset Item, itemset* result);
bitset<ITEM_SIZE> getItem(Tool T, Graph G, OwnStack S, Vertex V);
bool hasItem(itemset I, itemset Ilist);
void printGraph(OwnStack R, Tool T);
vector<int> adjList(BanList *Ban,  Graph G, OwnStack S);
void allPrint(OwnStack R);
double Pcmh(Param P, Tool T, Graph G, OwnStack S);
void writeSignificantsToFile(std::string filename, std::multimap<double, _OwnStack> container);
