#include "stat.hpp"
#include "igam.hpp"
#include "igami.hpp"

double MACHEP = 1.11022302462515654042E-16; // 2**-53
double MAXLOG = 8.8029691931113054295988E1; // log(2**127)

namespace stat {

// 0 <= p <= 1
double isf(double p) {
  double df = 1.0;
  double x = igamci(0.5 * df, p);
  return 2.0 * x;
}

double sf(double x) {
  double df = 1.0;
  return igamc(0.5 * df, 0.5 * x);
}

} // namespace stat