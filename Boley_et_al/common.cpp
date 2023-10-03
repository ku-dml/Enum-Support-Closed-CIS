/********** common.cpp **********/

/***** includes *****/
#include <cfloat>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <fstream>
#include <cmath>

// C++ libraries
#include <algorithm>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <vector>
using namespace std;

// custom libraries
#include "common.h"
#include "define.h"
#include "mylib.h"

/***** check the arguments *****/
void checkArgs(int argc, char *argv[]) {
  if (argc < 7) {
    fprintf(stderr,
            "usage: %s "
            "(pattern_file)(graph_file)(phenotype_file)(population_file)(theta)"
            "(alpha)[(param1)(value1)(param2)...]\n\n",
            argv[0]);
    fprintf(stderr, "list of params:\n");
    fprintf(stderr, "  -seed <INT> ... random seed (%d)\n", INI_seed);
    fprintf(stderr, "  -sigma <INT> ... minimum size of connector (%d)\n",
            INI_sigma);
    fprintf(stderr,
            "  -delta <DOUBLE> ... minimum density of a connector (%g)\n",
            INI_delta);
    fprintf(stderr,
            "  -tlim <DOUBLE> ... time limit in sec; if <0, it is inf. (%g)\n",
            INI_tlim);
    fprintf(stderr, "  -reduce <BOOL> ... whether reduction is invoked (%d)\n",
            INI_reduce);
    fprintf(stderr, "  -vtable <STR> ... file that describes relation between "
                    "vertex id & name; itable should be also specified\n");
    fprintf(stderr, "  -itable <STR> ... file that describes relation between "
                    "item id & name; vtable should be also specified\n");
    fprintf(
        stderr,
        "  -outname <STR> ... filename of output from vtable & itable (%s)\n",
        INI_outname);
    fprintf(stderr, "  -distname <STR> ... filename for distribution\n",
            INI_outname);
    fprintf(stderr, "  -ramub <INT> ... upper bound on used amount of RAM (%d)\n",
            INI_ramub);
    exit(EXIT_FAILURE);
  }
}

/***** read the arguments *****/
void readArgs(Param P, Graph G, int argc, char *argv[]) {
  char *arg;
  int k;
  P->ptn_file = argv[1];
  P->grh_file = argv[2];
  P->phenotype_file = argv[3];
  P->population_file = argv[4];
  P->theta = atoi(argv[5]);
  P->alpha = atof(argv[6]);
  P->turnWidth = INI_turnWidth;
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

  for (k = 7; k < argc; k += 2) {
    if (k == argc - 1 || argv[k][0] != '-') {
      fprintf(stderr, "error: arguments are not set appropriately.\n");
      exit(EXIT_FAILURE);
    }
    arg = argv[k] + 1;
    if (strcmp(arg, "seed") == Equiv)
      P->seed = atoi(argv[k + 1]);
    else if (strcmp(arg, "sigma") == Equiv)
      P->sigma = atoi(argv[k + 1]);
    else if (strcmp(arg, "delta") == Equiv)
      P->delta = atof(argv[k + 1]);
    else if (strcmp(arg, "turnWidth") == Equiv)
      P->turnWidth = atof(argv[k + 1]);
    else if (strcmp(arg, "tlim") == Equiv)
      P->tlim = atof(argv[k + 1]);
    else if (strcmp(arg, "reduce") == Equiv)
      P->reduce = (bool)atoi(argv[k + 1]);
    else if (strcmp(arg, "vtable") == Equiv)
      P->vtable = argv[k + 1];
    else if (strcmp(arg, "itable") == Equiv)
      P->itable = argv[k + 1];
    else if (strcmp(arg, "outname") == Equiv)
      strcpy(P->outname, argv[k + 1]);
    else if (strcmp(arg, "distname") == Equiv) {
      P->distname = new char[LEN_FILENAME];
      strcpy(P->distname, argv[k + 1]);
    } else if (strcmp(arg, "ramub") == Equiv)
      P->ramub = atoi(argv[k + 1]);
    else {
      fprintf(stderr, "error: <%s> is illegal parameter.\n", arg);
      exit(EXIT_FAILURE);
    }
  }
  if ((P->vtable == NULL && P->itable != NULL) ||
      (P->vtable != NULL && P->itable == NULL)) {
    fprintf(stderr, "error: if one of vtable and itable is specified, then the "
                    "other should be also specified.\n");
    exit(EXIT_FAILURE);
  }
}

/***** preprocess *****/
void preprocess(Param P, Tool T, Graph G) {
  char *str, *str_id, *str_item, **s;
  FILE *in;
  int b, i, j, items;

  str = new char[LINE_MAX];
  str_id = new char[LINE_MAX];
  str_item = new char[LINE_MAX];

  /* first scan of pattern file */
  G->n = 0;
  // G->m = 0;
  G->items = 0;
  G->item_density = 0.;
  in = open_file(P->ptn_file, "r");
  while (fgets(str, LINE_MAX - 1, in) != NULL) {
    if (str[0] == '#')
      continue;
    sscanf(str, "%s%s", str_id, str_item);
    auto itr = T->VMap.find(atoi(str_id));
    if (itr == T->VMap.end()) { // if str_id is not in the VMap
      T->VMap[atoi(str_id)] = G->n;
      T->VMapInv[G->n] = atoi(str_id);
      G->n++;
    }
    s = split(str_item, ',');
    items = getNumChar(str_item, ',') + 1;
    for (j = 0; j < items; j++) {
      i = atoi(s[j]);
      auto itr = T->IMap.find(i);
      if (itr == T->IMap.end()) { // if i is not in the IMap
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

  T->reduced_vertices = 0;
  T->reduced_edges = 0;
  for (j = 0; j < G->items; j++) {
    VertexIDSeq empty;
    T->seq_item.push_back(empty);
  }
  in = open_file(P->population_file, "r");
  T->numPop = 0;
  while (fgets(str, LINE_MAX - 1, in) != NULL) {
    if (str[0] == '#')
      continue;
    sscanf(str, "%s%s", str_item, str_id);
    auto itr = T->PMap.find(atoi(str_id));
    if (itr == T->PMap.end()) { // if str_id is not in the PMap
      T->PMap[atoi(str_id)] = T->numPop;
      T->PMapInv[T->numPop] = atoi(str_id);
      T->numPop++;
    }
  }
  fclose(in);
  itemset I;
  I.reset();
  for (int i = 0; i < T->numPop; i++) {
    T->population.push_back(I);
  }
  delete[] str;
  delete[] str_id;
  delete[] str_item;
}

void initTool(Tool T, Graph G) {
  int i, j;
  for (i = 0; i < G->n; i++)
    if (G->V[i] != NULL)
      T->seq_all.push_back(i);
}

/***** read the pattern file *****/
void readPattern(Param P, Tool T, Graph G) {
  char *str, *str_id, *str_item, **s;
  FILE *in;
  double item_density = 0.;
  int x, i, j, p;

  str = new char[LINE_MAX];
  str_id = new char[LINE_MAX];
  str_item = new char[LINE_MAX];
  in = open_file(P->ptn_file, "r");

  while (fgets(str, LINE_MAX - 1, in) != NULL) {
    if (str[0] == '#')
      // comment?
      continue;
    strcpy(str_item, ""); // for the case itemset is empty
    sscanf(str, "%s%s", str_id, str_item);
    x = T->VMap[atoi(str_id)];
    G->V[x] = new struct _Vertex;
    G->V[x]->id = x;
    G->V[x]->deg = 0;
    G->V[x]->Itv_ptr = new IntvIDSeq[1];
    G->V[x]->I = new itemset;
    G->V[x]->I->reset();

    if (strcmp(str_item, "") == Equiv || strcmp(str_item, "!") == Equiv) {
      G->V[x]->items = 0;
    } else {
      s = split(str_item, ',');
      G->V[x]->items = getNumChar(str_item, ',') + 1;
      item_density += (double)(G->V[x]->items);
      for (j = 0; j < G->V[x]->items; j++) {
        i = T->IMap[atoi(s[j])];
        G->V[x]->I->set(i);
        T->seq_item[i].push_back(x);
        free(s[j]);
      }
      free(s);
    }
    /* Reduction 1 */
    if (P->reduce && G->V[x]->items < P->theta) {
      delete G->V[x]->I;
      delete G->V[x];
      G->V[x] = NULL;
      T->reduced_vertices++;
    }
  }

  G->item_density = (item_density / (double)G->items) / (double)G->n;

  delete[] str;
  delete[] str_id;
  delete[] str_item;
  fclose(in);
}

void readPhenotype(Param P, Tool T) {
  char *str, *str_iid, *str_phenotype;
  FILE *in;
  int i = 0;
  str = new char[LINE_MAX];
  str_iid = new char[LINE_MAX];
  str_phenotype = new char[LINE_MAX];
  in = open_file(P->phenotype_file, "r");
  T->phenotype.reset();
  while (fgets(str, LINE_MAX - 1, in) != NULL) {
    if (str[0] == 'I') {
      // skip first line
      continue;
    }
    sscanf(str, "#%s%s", str_iid, str_phenotype);
    i = T->IMap[atoi(str_iid) - 1];
    if (atoi(str_phenotype) == 1) {
      T->phenotype.set(i);
    }
    i++;
  }
  delete[] str;
  delete[] str_iid;
  delete[] str_phenotype;
  fclose(in);
}

void readPopulation(Param P, Tool T) {
  char *str, *str_iid, *str_population;
  FILE *in;
  int i, pop;
  str = new char[LINE_MAX];
  str_iid = new char[LINE_MAX];
  str_population = new char[LINE_MAX];
  in = open_file(P->population_file, "r");
  itemset I;
  I.reset();
  while (fgets(str, LINE_MAX - 1, in) != NULL) {
    if (str[0] == '#')
      continue;
    sscanf(str, "%s%s", str_iid, str_population);
    i = T->IMap[atoi(str_iid)-1];
    pop = T->PMap[atoi(str_population)];
    T->population[pop].set(i);
  }
  delete[] str;
  delete[] str_iid;
  delete[] str_population;
  fclose(in);
}

void computeStatisticVal(Tool T) {
  for (int i = 0; i < T->numPop; i++) {
    T->n.push_back(T->population[i].count());
    T->n1.push_back((T->phenotype & T->population[i]).count());
    T->n2.push_back(T->n[i] - T->n1[i]);
  }
}

/***** read the graph file *****/
void readGraph(Param P, Tool T, Graph G) {
  char *str;
  FILE *in;
  itemset I;
  double w;
  int a, b, x, y, line = 0;

  str = new char[LINE_MAX];
  in = open_file(P->grh_file, "r");

  while (fgets(str, LINE_MAX - 1, in) != NULL) {
    line++;
    if (str[0] == '#') // comment is skipped
      continue;
    sscanf(str, "%d%lf%d", &a, &w, &b);
    if (a == b) {
      fprintf(stderr, "The line %d is ignored since %d=%d (self-loop).\n", line,
              a, b);
      continue;
    }
    x = T->VMap[a];
    y = T->VMap[b];
    if (P->reduce) {
      if (G->V[x] == NULL || G->V[y] == NULL)
        continue;
      /* reduction 2 */
      int p, items = 0;
      I = *(G->V[x]->I) & *(G->V[y]->I);
      items = I.count();
      if (items < P->theta) {
        T->reduced_edges++;
        continue;
      }
    }

    bool add = true;
    for (auto itr = G->V[x]->A.begin(); itr != G->V[x]->A.end(); itr++) {
      if (*itr == y) {
        add = false;
        break;
      }
    }
    if (add)
      G->V[x]->A.push_back(y);
  }
  G->m = countNumEdges(G);
  delete[] str;
  fclose(in);
}

/***** count the number of edges *****/
int countNumEdges(Graph G) {
  int m = 0;
  for (int x = 0; x < G->n; x++) {
    if (G->V[x] == NULL)
      continue;
    int deg = 0;
    for (auto adj = G->V[x]->A.begin(); adj != G->V[x]->A.end(); adj++)
      deg++;
    m += deg;
  }
  if (m % 2 == 1)
    fprintf(stderr, "warning: degree sum is not even.\n");
  return m / 2;
}

/***** prepare BFS *****/
unsigned int prepareBFS(BFSTool B, int n, VertexIDSeq &seq) {
  unsigned int marker;
  int size;
  marker = B->marker;
  B->marker++;
  // reset the marks of all vertices if marker=marker_max
  if (B->marker > B->marker_max - 1) {
    for (int x = 0; x < n; x++)
      B->mark[x] = 0;
    B->marker = 1;
  }
  size = seq.size();
  for (int i = 0; i < size; i++) {
    if (seq[i] != -1) {
      B->mark[seq[i]] = marker;
    }
  }
  return marker;
}

int getClosure(Param P, Tool T, Graph G, BFSTool B, Solution *S, BanList *Ban,
               Vertex V, itemset Item, itemset *result) {
  itemset I = Item & *(V->I);
  int numItems = I.count();
  if (numItems < P->theta)
    return 0;
  else if (Ban->bit[V->id])
    return 0;
  queue<int> Q;
  int top, nextV;
  Q.push(V->id);
  auto end_itr = S->seq.end();
  for (auto itr = S->seq.begin(); itr != end_itr; ++itr) {
    if (*itr != -1)
      Q.push(G->V[*itr]->id);
  }
  unsigned int marker = prepareBFS(B, G->n, S->seq);

  B->marker =
      B->marker + 1; // mark of v which was push on Q but not checked about C
  while (!Q.empty()) {
    top = Q.front();
    if (B->mark[top] != marker) {
      S->push_back(top);
      B->mark[top] = marker;
    }
    auto end = G->V[top]->A.end();
    for (auto itr = G->V[top]->A.begin(); itr != end; itr++) {
      nextV = *itr;
      if (B->mark[nextV] < marker && hasItem(*(G->V[nextV]->I), I)) {
        if (Ban->bit[nextV])
          return 0;
        Q.push(nextV);
        B->mark[nextV] = marker + 1;
      }
    }
    Q.pop();
  }
  *result = I;
  if (numItems == P->theta)
    return 2;
  return 1;
}

itemset getItem(Tool T, Graph G, OwnStack S, Vertex V) {
  itemset I;
  I.set();
  for (auto v_itr = S->seq.begin(); v_itr != S->seq.end(); v_itr++) {
    if (*v_itr != -1)
      I = I & *(G->V[*v_itr]->I);
    if (I.none())
      break;
  }
  I = I & *(V->I);
  return I;
}

bool hasItem(itemset I, itemset Ilist) { return ((I & Ilist) == Ilist); }

void printGraph(OwnStack C, Tool T) {
  auto end = C->seq.end();
  for (auto itr = C->seq.begin(); itr != end; ++itr) {
    if (*itr != -1) {
      cout << T->VMapInv[*itr] << " ";
    }
  }
}

vector<int> adjList(BanList *Ban, Graph G, OwnStack S) {
  vector<int> adjVList;
  vertexset check;
  auto end = S->seq.end();
  int id;
  for (auto itr = S->seq.begin(); itr != end; ++itr) {
    if (*itr == -1)
      continue;
    for (auto itr2 = G->V[*itr]->A.begin(); itr2 != G->V[*itr]->A.end();
         itr2++) {
      id = *itr2;
      if (check[id] || S->bit[id] || Ban->bit[id])
        continue;
      adjVList.push_back(id);
      check.set(id);
    }
  }
  return adjVList;
}

double Pcmh(Param P, Tool T, Graph G, OwnStack S) {
  itemset I, sepI;
  int nj, n1j, n2j, aSj, xSj;
  double denom = 0.0, numer = 0.0;
  I.set();
  for (auto v_itr = S->seq.begin(); v_itr != S->seq.end(); v_itr++) {
    if (*v_itr != -1)
      I = I & *(G->V[*v_itr]->I);
    if (I.none())
      break;
  }
  for (int i = 0; i < T->numPop; i++) {
    xSj = ((~I) & T->population[i]).count();
    aSj = ((~I) & T->population[i] & T->phenotype).count();
    numer += 1.0 * aSj - (xSj * (1.0 * T->n1[i] / T->n[i]));
    denom += (1.0 * T->n1[i] * T->n2[i] * xSj * (T->n[i] - xSj)) /
             ((T->n[i] - 1.0) * T->n[i] * T->n[i]);
  }

  numer = numer * numer;
  return 1.0 - erf(sqrt(0.5 * numer / denom));
}

void writeSignificantsToFile(std::string filename, std::multimap<double, _OwnStack> container) {
  cout << "writing to file... ";

  std::ofstream file;
  file.open(filename, std::ios::out);
  file << "p-value,v-indices" << std::endl;
  auto end_itr = container.end();
  for (auto c_itr = container.begin(); c_itr != end_itr; ++c_itr) {
    // write p-value,
    file << c_itr->first;
    // then, v
    for (auto &e: c_itr->second.seq) {
      file << "," << e;
    }
    file << std::endl;
  }


 file.close();
  cout << "finished!" << endl;
}
