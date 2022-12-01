/*
 * Cephes Math Library Release 2.1:  December, 1988
 * Copyright 1984, 1987, 1988 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

/* Sources:
 * [1] Holin et. al., "Polynomial and Rational Function Evaluation",
 *     https://www.boost.org/doc/libs/1_61_0/libs/math/doc/html/math_toolkit/roots/rational.html
 */

/* Scipy changes:
 * - 06-23-2016: add code for evaluating rational functions
 */

#include "polevl.hpp"
#include <cmath>

// static NPY_INLINE double polevl(double x, const double coef[], int N)
double polevl(double x, const double coef[], int N_pol) {
  double ans;
  int i;
  const double *p;

  p = coef;
  ans = *p++;
  i = N_pol;

  do {
    ans = ans * x + *p++;
  } while (--i);

  return ans;
}

/*                                                     p1evl() */
/*                                          N
 * Evaluate polynomial when coefficient of x  is 1.0.
 * Otherwise same as polevl.
 */

// static NPY_INLINE double p1evl(double x, const double coef[], int N)
double p1evl(double x, const double coef[], int N_pol) {
  double ans;
  const double *p;
  int i;

  p = coef;
  ans = x + *p++;
  i = N_pol - 1;

  do
    ans = ans * x + *p++;
  while (--i);

  return (ans);
}

/* Evaluate a rational function. See [1]. */

// static NPY_INLINE double ratevl(double x, const double num[], int M,
double ratevl(double x, const double num[], int M,
                            const double denom[], int N_pol) {
  int i, dir;
  double y, num_ans, denom_ans;
  double absx = fabs(x);
  const double *p;

  if (absx > 1) {
    /* Evaluate as a polynomial in 1/x. */
    dir = -1;
    p = num + M;
    y = 1 / x;
  } else {
    dir = 1;
    p = num;
    y = x;
  }

  /* Evaluate the numerator */
  num_ans = *p;
  p += dir;
  for (i = 1; i <= M; i++) {
    num_ans = num_ans * y + *p;
    p += dir;
  }

  /* Evaluate the denominator */
  if (absx > 1) {
    p = denom + N_pol;
  } else {
    p = denom;
  }

  denom_ans = *p;
  p += dir;
  for (i = 1; i <= N_pol; i++) {
    denom_ans = denom_ans * y + *p;
    p += dir;
  }

  if (absx > 1) {
    i = N_pol - M;
    return pow(x, i) * num_ans / denom_ans;
  } else {
    return num_ans / denom_ans;
  }
}
