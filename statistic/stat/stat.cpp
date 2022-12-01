#include "stat.hpp"
#include "igami.hpp"

double MACHEP = 1.11022302462515654042E-16; // 2**-53
double MAXLOG = 8.8029691931113054295988E1; // log(2**127)

// 0 <= p <= 1
double isf(double p) {
  double df = 1.0;
  double x = igamci(0.5 * df, p);
  return 2.0 * x;
}