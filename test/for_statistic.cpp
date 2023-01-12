#include "../statistic/p_value.hpp"
#include "../statistic/stat/stat.hpp"
#include <assert.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <tuple>
#include <vector>

double ERR = 1e-10;

void assert_eq(double a, double b);
void test1();
void test2();
void test3();
void test4();
void test5();
double for_minimal_p_value(std::vector<double> xs, std::vector<double> n1s,
                           std::vector<double> n2s);
double for_envelope(std::vector<double> xs, std::vector<double> n1s,
                    std::vector<double> n2s);
double for_p_value(std::vector<double> xs, std::vector<double> n1s,
                   std::vector<double> n2s, std::vector<double> as);

_Tool testcase1_tool();
std::tuple<itemset, itemset, itemset, itemset, itemset> testcase1_items();

int main() {
  test1();
  test2();
  test3();
  test4();
  test5();

  printf("\nAll test passed!\n");
  return 0;
}

void assert_eq(double a, double b) {
  printf("test: %f == %f\n", a, b);
  assert(fabs(a - b) < ERR);
}

// for inverse_threshold
void test1() {
  printf("\nTest 1:\n");
  _Tool T;
  T.numPop = 1;
  _Stat stat = _Stat(&T, nullptr);
  std::vector<double> testcase1{
      0.7,
      0.4401044134621469,
      0.6914818268791262,
      0.9548036573918648,
      0.04222840107181414,
      0.3055911067795858,
      0.234451381603198,
      0.024522097396397324,
      0.4190335107835975,
      0.49135687067728073,
      0.8574504652721333,
      0.25891434778568556,
      0.8022280733954448,
      0.2857279096602334,
      0.26451967505298934,
  };
  std::vector<double> expected1{
      0.14847186183254535,  0.5960101059959402,  0.1574858868158796,
      0.003212117054733589, 4.126026230109639,   1.049635250994285,
      1.413650209588385,    5.057322338386889,   0.6530242253250768,
      0.473553604963164,    0.03226354437758338, 1.2745526092845478,
      0.0627325185109155,   1.139642367880397,   1.2449475050363494,
  };
  for (int i = 0; i < testcase1.size(); ++i) {
    assert_eq(stat.inverse_threshold(testcase1[i], 1), expected1[i]);
  }
}

// for minimal_p_value
void test2() {
  printf("\nTest 2:\n");
  // case 1
  _Tool T = testcase1_tool();
  _Graph G;

  _Vertex *V = (_Vertex *)malloc(5 * sizeof(_Vertex));
  auto items = testcase1_items();

  _Vertex v1, v2, v3, v4, v5;
  v1.I = &std::get<0>(items);
  V[0] = v1;
  v2.I = &std::get<1>(items);
  V[1] = v2;
  v3.I = &std::get<2>(items);
  V[2] = v3;
  v4.I = &std::get<3>(items);
  V[3] = v4;
  v5.I = &std::get<4>(items);
  V[4] = v5;

  G.V = &V;
  _Stat stat(&T, &G);
  _OwnStack S;
  std::vector<int> v{0, 1, 2};
  S.seq = v;
  std::vector<double> xs1{2.0, 1.0, 3.0};
  std::vector<double> n1s1{2.0, 1.0, 3.0};
  std::vector<double> n2s1{1.0, 1.0, 1.0};
  assert_eq(for_minimal_p_value(xs1, n1s1, n2s1), 5.568421052631581);
  assert_eq(stat.minimal_p_value(&S), 5.568421052631581);

  // case 2
  std::vector<double> xs2{2.0, 2.0, 2.0, 2.0};
  std::vector<double> n1s2{1.0, 1.0, 1.0, 2.0};
  std::vector<double> n2s2{1.0, 1.0, 1.0, 0.0};
  assert_eq(for_minimal_p_value(xs2, n1s2, n2s2), 0.0);

  // case 3
  std::vector<double> xs3{2.0, 3.0, 4.0, 5.0, 2.0, 3.0, 4.0};
  std::vector<double> n1s3{3.0, 1.0, 2.0, 3.0, 1.0, 0.0, 1.0};
  std::vector<double> n2s3{3.0, 2.0, 5.0, 6.0, 2.0, 4.0, 3.0};
  assert_eq(for_minimal_p_value(xs3, n1s3, n2s3), 12.633686016585644);
}

// for envelope
void test3() {
  printf("\nTest 3:\n");
  // case 1
  _Tool T = testcase1_tool();
  _Graph G;

  _Vertex *V = (_Vertex *)malloc(5 * sizeof(_Vertex));
  auto items = testcase1_items();

  _Vertex v1, v2, v3, v4, v5;
  v1.I = &std::get<0>(items);
  V[0] = v1;
  v2.I = &std::get<1>(items);
  V[1] = v2;
  v3.I = &std::get<2>(items);
  V[2] = v3;
  v4.I = &std::get<3>(items);
  V[3] = v4;
  v5.I = &std::get<4>(items);
  V[4] = v5;

  G.V = &V;
  _Stat stat(&T, &G);
  _OwnStack S;
  std::vector<int> v{0, 1, 2};
  S.seq = v;
  std::vector<double> xs1{2.0, 1.0, 3.0};
  std::vector<double> n1s1{2.0, 1.0, 3.0};
  std::vector<double> n2s1{1.0, 1.0, 1.0};
  assert_eq(for_envelope(xs1, n1s1, n2s1), 5.568421052631581);
  assert_eq(stat.envelope(&S), 5.568421052631581);

  // case 2
  std::vector<double> xs2{2.0, 2.0, 2.0, 2.0};
  std::vector<double> n1s2{1.0, 1.0, 1.0, 2.0};
  std::vector<double> n2s2{1.0, 1.0, 1.0, 0.0};
  assert_eq(for_envelope(xs2, n1s2, n2s2), 0.0);

  // case 3
  std::vector<double> xs3{2.0, 3.0, 4.0, 5.0, 2.0, 3.0, 4.0};
  std::vector<double> n1s3{3.0, 1.0, 2.0, 3.0, 1.0, 0.0, 1.0};
  std::vector<double> n2s3{3.0, 2.0, 5.0, 6.0, 2.0, 4.0, 3.0};
  assert_eq(for_envelope(xs3, n1s3, n2s3), 20.70057725466677);
}

// for p_value
void test4() {
  printf("\nTest 4:\n");
  // case 1
  _Tool T = testcase1_tool();
  _Graph G;

  _Vertex *V = (_Vertex *)malloc(5 * sizeof(_Vertex));
  auto items = testcase1_items();

  _Vertex v1, v2, v3, v4, v5;
  v1.I = &std::get<0>(items);
  V[0] = v1;
  v2.I = &std::get<1>(items);
  V[1] = v2;
  v3.I = &std::get<2>(items);
  V[2] = v3;
  v4.I = &std::get<3>(items);
  V[3] = v4;
  v5.I = &std::get<4>(items);
  V[4] = v5;

  G.V = &V;
  _Stat stat(&T, &G);
  _OwnStack S;
  std::vector<int> v{0, 1, 2};
  S.seq = v;
  std::vector<double> xs1{2.0, 1.0, 3.0};
  std::vector<double> n1s1{2.0, 1.0, 3.0};
  std::vector<double> n2s1{1.0, 1.0, 1.0};
  std::vector<double> as{1.0, 1.0, 2.0};
  assert_eq(for_p_value(xs1, n1s1, n2s1, as), 0.01052631578947361);
  assert_eq(stat.p_value(&S), 0.01052631578947361);

  // case 2
  std::vector<double> xs2{2.0, 2.0, 2.0, 2.0};
  std::vector<double> n1s2{1.0, 1.0, 1.0, 2.0};
  std::vector<double> n2s2{1.0, 1.0, 1.0, 0.0};
  assert_eq(for_p_value(xs2, n1s2, n2s2, as), 0.0);

  // case 3
  std::vector<double> xs3{2.0, 3.0, 4.0, 5.0, 2.0, 3.0, 4.0};
  std::vector<double> n1s3{3.0, 1.0, 2.0, 3.0, 1.0, 0.0, 1.0};
  std::vector<double> n2s3{3.0, 2.0, 5.0, 6.0, 2.0, 4.0, 3.0};
  std::vector<double> as2{1.0, 2.0, 4.0, 3.0, 1.0, 1.0, 2.0};
  assert_eq(for_p_value(xs3, n1s3, n2s3, as2), 35.69345152988275);
}

double for_minimal_p_value(std::vector<double> xs, std::vector<double> n1s,
                           std::vector<double> n2s) {
  double num = 0.0, denum = 0.0;
  double n1sj, n2sj, xsj, nsj;
  double aMin = 0.0, aMax = 0.0;
  for (int j = 0; j < xs.size(); ++j) {
    xsj = xs[j];
    n1sj = n1s[j];
    n2sj = n2s[j];
    nsj = n1sj + n2sj;

    num += xsj * n1sj / nsj;
    denum += n1sj * n2sj * xsj * (1.0 - xsj / nsj) / (nsj * (nsj - 1.0));

    aMin += std::max(0.0, xsj - n2sj);
    aMax += std::min(xsj, n1sj);
  }

  if (denum == 0.0) {
    return 0.0;
  } else {
    double numMin = (aMin - num) * (aMin - num);
    double numMax = (aMax - num) * (aMax - num);
    return std::max(numMin, numMax) / denum;
  }
}

double for_envelope(std::vector<double> xs, std::vector<double> n1s,
                    std::vector<double> n2s) {
  int J = xs.size();

  double n1sj, n2sj, xsj, nsj, nsj_2;
  std::vector<std::tuple<double, int>> betaLs(J), betaRs(J);
  std::vector<double> xStarsL(J), xStarsR(J);
  for (int j = 0; j < J; ++j) {
    xsj = xs[j];
    n1sj = (double)n1s[j];
    n2sj = (double)n2s[j];
    nsj = n1sj + n2sj;
    nsj_2 = nsj * nsj;
    xStarsL[j] = nsj;
    xStarsR[j] = nsj;

    betaLs[j] = std::make_tuple(n1sj * std::max(xsj, n2sj) / nsj_2, j);
    betaRs[j] = std::make_tuple(n2sj * std::max(xsj, n1sj) / nsj_2, j);
  }
  std::sort(betaLs.begin(), betaLs.end());
  std::sort(betaRs.begin(), betaRs.end());

  double env = 0.0;
  int piL, piR;
  for (int j = 0; j < J; ++j) {
    piL = std::get<1>(betaLs[j]);
    piR = std::get<1>(betaRs[j]);

    xStarsL[piL] = std::max(xs[piL], (double)n2s[piL]);
    env = std::max(env, for_minimal_p_value(xStarsL, n1s, n2s));

    xStarsR[piR] = std::max(xs[piL], (double)n1s[piL]);
    env = std::max(env, for_minimal_p_value(xStarsR, n1s, n2s));
  }

  return env;
}

double for_p_value(std::vector<double> xs, std::vector<double> n1s,
                   std::vector<double> n2s, std::vector<double> as) {
  double a = 0.0, num = 0.0, denum = 0.0;
  double n1sj, n2sj, xsj, nsj, asj;
  for (int j = 0; j < xs.size(); ++j) {
    xsj = (double)xs[j];
    asj = (double)as[j];
    n1sj = (double)n1s[j];
    n2sj = (double)n2s[j];
    nsj = (double)n1sj + n2sj;

    a += asj;
    num += xsj * n1sj / nsj;
    denum += n1sj * n2sj * xsj * (1.0 - xsj / nsj) / (nsj * (nsj - 1.0));
  }

  if (denum == 0.0) {
    return 0.0;
  } else {
    num = a - num;
    num *= num;
    return num / denum;
  }
}

std::tuple<itemset, itemset, itemset, itemset, itemset> testcase1_items() {
  itemset i1;
  i1.set();
  i1.reset(0);
  i1.reset(2);
  i1.reset(3);

  itemset i2;
  i2.set();
  i2.reset(5);

  itemset i3;
  i3.set();
  i3.reset(6);
  i3.reset(8);

  itemset i4;
  i4.set();
  i4.reset(7);

  itemset i5;
  i5.set();
  i5.reset(1);
  i5.reset(4);

  return std::make_tuple(i1, i2, i3, i4, i5);
}

_Tool testcase1_tool() {
  _Tool T;

  T.numPop = 3;

  std::vector<itemset> population(T.numPop);

  itemset i1;
  i1.reset();
  i1.set(0);
  i1.set(1);
  i1.set(2);
  population[0] = i1;

  itemset i2;
  i2.reset();
  i2.set(3);
  i2.set(4);
  population[1] = i2;

  itemset i3;
  i3.reset();
  i3.set(5);
  i3.set(6);
  i3.set(7);
  i3.set(8);
  population[2] = i3;

  T.population = population;

  itemset phenotype;
  phenotype.reset();
  phenotype.set(0);
  phenotype.set(1);
  phenotype.set(3);
  phenotype.set(4);
  phenotype.set(5);
  phenotype.set(6);
  phenotype.set(7);

  T.phenotype = phenotype;

  std::vector<int> n{3, 2, 4};
  std::vector<int> n1{2, 1, 3};
  std::vector<int> n2{1, 1, 1};

  T.n = n;
  T.n1 = n1;
  T.n2 = n2;

  return T;
}

// for inverse_threshold
void test5() {
  printf("\nTest 5:\n");
  _Tool T;
  T.numPop = 1;
  _Stat stat = _Stat(&T, nullptr);
  std::vector<double> testcase5{0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
                                0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0};
  std::vector<double> expected5{
      0.47950012218695337, 0.7518296340458492,  0.6547208460185768,
      0.583882420770365,   0.5270892568655381,  0.47950012218695337,
      0.4385780260809997,  0.40278369424647564, 0.37109336952269756,
      0.34278171114790873, 0.31731050786291115, 0.15729920705028105,
      0.08326451666355042, 0.04550026389635857,
  };
  for (int i = 0; i < testcase5.size(); ++i) {
    assert_eq(stat.survival_function(testcase5[i]), expected5[i]);
  }
}
