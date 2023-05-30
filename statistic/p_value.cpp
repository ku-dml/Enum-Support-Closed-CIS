#include "p_value.hpp"
#include "../Boley_et_al/define.h"
#include "stat/stat.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

_Stat::_Stat(Tool T, Graph G) {
  // copy address
  this->tool = T;
  this->graph = G;

  int J = T->numPop;
  this->J = J;
  // reserve buffer
  itemset I;
  this->I_buffer = std::move(I);
  this->xs_buffer = std::vector<double>(J);
  this->betaLs_buffer = std::vector<std::tuple<double, int>>(J);
  this->betaRs_buffer = std::vector<std::tuple<double, int>>(J);
  this->xStarsL_buffer = std::vector<double>(J);
  this->xStarsR_buffer = std::vector<double>(J);
}

double _Stat::inverse_threshold(double alpha, int k) {
  return stat::isf(alpha / (double)k);
}

double _Stat::survival_function(double x) { return stat::sf(x); }

double _Stat::p_value(OwnStack S) {
  // construct item set in order to compute xs
  this->I_buffer.set();
  auto end_itr = S->seq.end();
  for (auto v_itr = S->seq.begin(); v_itr != end_itr; ++v_itr) {
    if (*v_itr != -1) {
      this->I_buffer &= *(this->graph->V[*v_itr]->I);
    }
    if (this->I_buffer.none()) {
      break;
    }
  }
  // cout << ~this->I_buffer << endl;
  // cout << (~this->I_buffer).count() << endl;
  // auto v = std::vector<int>();
  // auto a_tmp = 0;
  // for (int i = 0; i < this->I_buffer.size(); ++i) {
  //   if (this->I_buffer[i] == 1) {
  //     v.push_back(this->tool->IMapInv[i]);
  //     if (this->tool->phenotype[i] == 1) {
  //       a_tmp++;
  //     }
  //     // printf("%d, ", this->tool->IMapInv[i]);
  //   }
  // }
  // std::sort(v.begin(), v.end());
  // for (auto &e: v) {
  //   cout << e << ", ";
  // }
  // cout << " -> " << a_tmp << endl;


  double a = 0.0, num = 0.0, denum = 0.0;
  double n1sj, n2sj, xsj, nsj, asj;
  for (int j = 0; j < this->J; ++j) {
    xsj = (double)((~(this->I_buffer)) & this->tool->population[j]).count();
    asj = (double)(((~(this->I_buffer)) & this->tool->population[j]) &
                   this->tool->phenotype)
              .count();
    n1sj = (double)(this->tool->n1[j]);
    n2sj = (double)(this->tool->n2[j]);
    nsj = (double)(this->tool->n[j]);

    // cout << asj << " " << xsj << " " << n1sj << " " << n2sj << endl;
    // cout << ~this->I_buffer << endl;
    // cout << this->tool->population[0] << endl;
    // cout << this->tool->phenotype << endl;

    // printf("\n%f, %f, %f, %f, %f\n", xsj, asj, n1sj, n2sj, nsj);
    // cout << ((~this->I_buffer) & this->tool->population[j] & this->tool->phenotype).count() << endl;

    a += asj;
    num += xsj * n1sj / nsj;
    denum += n1sj * n2sj * xsj * (1.0 - xsj / nsj) / (nsj * (nsj - 1.0));
  }

  if (denum == 0.0) {
    // return 0.0;
    // return __DBL_MAX__;
  } else {
    num = a - num;
    num *= num;
    return num / denum;
  }
}

double _Stat::minimal_p_value(OwnStack S) {
  // construct item set in order to compute xs
  this->I_buffer.set();

  auto end_itr = S->seq.end();
  for (auto v_itr = S->seq.begin(); v_itr != end_itr; ++v_itr) {
    if (*v_itr != -1) {
      this->I_buffer &= *(this->graph->V[*v_itr]->I);
    }
    if (this->I_buffer.none()) {
      break;
    }
  }

  double num = 0.0, denum = 0.0;
  double n1sj, n2sj, xsj, nsj;
  double aMin = 0.0, aMax = 0.0;
  for (int j = 0; j < this->J; ++j) {
    xsj = (double)((~(this->I_buffer)) & this->tool->population[j]).count();
    n1sj = (double)(this->tool->n1[j]);
    n2sj = (double)(this->tool->n2[j]);
    nsj = (double)(this->tool->n[j]);


    num += xsj * n1sj / nsj;
    denum += n1sj * n2sj * xsj * (1.0 - xsj / nsj) / (nsj * (nsj - 1.0));

    aMin += std::max(0.0, xsj - n2sj);
    aMax += std::min(xsj, n1sj);
  }
  if (denum == 0.0) {
    return 0.0;
  } else {
    double numMin = aMin - num;
    numMin *= numMin;
    double numMax = aMax - num;
    numMax *= numMax;
    return std::max(numMin, numMax) / denum;
  }
}

// This private function is only used in envelope function
// in order to compute minimal p-value from array.
double _Stat::minimal_p_value_inner(const std::vector<double> &xs) {
  double num = 0.0, denum = 0.0;
  double n1sj, n2sj, xsj, nsj;
  double aMin = 0.0, aMax = 0.0;
  for (int j = 0; j < this->J; ++j) {
    xsj = xs[j];
    n1sj = (double)(this->tool->n1[j]);
    n2sj = (double)(this->tool->n2[j]);
    nsj = n1sj + n2sj;

    // cout << xsj << " " << n1sj << " " << n2sj << " " << nsj << endl; 

    num += xsj * n1sj / nsj;
    denum += n1sj * n2sj * xsj * (1.0 - xsj / nsj) / (nsj * (nsj - 1.0));

    aMin += std::max(0.0, xsj - n2sj);
    aMax += std::min(xsj, n1sj);
  }

  if (denum == 0.0) {
    return 0.0;
  } else {
    double numMin = aMin - num;
    numMin *= numMin;
    double numMax = aMax - num;
    numMax *= numMax;
    return std::max(numMin, numMax) / denum;
  }
}

double _Stat::envelope(OwnStack S) {
  // cout << "endvelope" << endl;
  // cout << this->betaLs_buffer.size() << endl;
  // cout << this->J << endl;

  // construct item set in order to compute xs
  this->I_buffer.set();
  auto end_itr = S->seq.end();
  for (auto v_itr = S->seq.begin(); v_itr != end_itr; ++v_itr) {
    if (*v_itr != -1) {
      this->I_buffer &= *(this->graph->V[*v_itr]->I);
    }
    if (this->I_buffer.none()) {
      break;
    }
  }

  // then compute xs (use xs later)
  double n1sj, n2sj, xsj, nsj, nsj_2;
  for (int j = 0; j < this->J; ++j) {
    xsj = (double)((~(this->I_buffer)) & this->tool->population[j]).count();
    this->xs_buffer[j] = xsj;
    n1sj = (double)(this->tool->n1[j]);
    n2sj = (double)(this->tool->n2[j]);
    nsj = (double)(this->tool->n[j]);
    nsj_2 = nsj * nsj;
    this->xStarsL_buffer[j] = nsj;
    this->xStarsR_buffer[j] = nsj;

    this->betaLs_buffer[j] =
        std::make_tuple(n1sj * std::max(xsj, n2sj) / nsj_2, j);
    this->betaRs_buffer[j] =
        std::make_tuple(n2sj * std::max(xsj, n1sj) / nsj_2, j);
  }
  std::sort(this->betaLs_buffer.begin(), this->betaLs_buffer.end());
  std::sort(this->betaRs_buffer.begin(), this->betaRs_buffer.end());

  double env = 0.0;
  int piL, piR;
  for (int j = 0; j < this->J; ++j) {
    // cout << "for" << endl;
    piL = std::get<1>(this->betaLs_buffer[j]);
    piR = std::get<1>(this->betaRs_buffer[j]);

    this->xStarsL_buffer[piL] =
        std::max(this->xs_buffer[piL], (double)(this->tool->n2[piL]));
    env = std::max(env, this->minimal_p_value_inner(this->xStarsL_buffer));
    // cout << "1: " << this->minimal_p_value_inner(this->xStarsL_buffer) << endl;
    // for (auto e: this->xStarsL_buffer) {
    //   cout << e << endl;
    // }

    this->xStarsR_buffer[piR] =
        std::max(this->xs_buffer[piL], (double)(this->tool->n1[piL]));
    env = std::max(env, this->minimal_p_value_inner(this->xStarsR_buffer));
    // cout << "2: " << this->minimal_p_value_inner(this->xStarsR_buffer) << endl;
    // for (auto e: this->xStarsR_buffer) {
    //   cout << e << endl;
    // }
  }

  return env;
}
