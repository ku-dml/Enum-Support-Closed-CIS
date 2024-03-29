/*
 * Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
All rights reserved.
 */
extern double MACHEP;

#include "zeta.hpp"
#include <cmath>
#include <stdio.h>

/* Expansion coefficients
 * for Euler-Maclaurin summation formula
 * (2k)! / B2k
 * where B2k are Bernoulli numbers
 */
static double A_ZETA[] = {
    12.0,
    -720.0,
    30240.0,
    -1209600.0,
    47900160.0,
    -1.8924375803183791606e9, /*1.307674368e12/691 */
    7.47242496e10,
    -2.950130727918164224e12,  /*1.067062284288e16/3617 */
    1.1646782814350067249e14,  /*5.109094217170944e18/43867 */
    -4.5979787224074726105e15, /*8.028576626982912e20/174611 */
    1.8152105401943546773e17,  /*1.5511210043330985984e23/854513 */
    -7.1661652561756670113e18  /*1.6938241367317436694528e27/236364091 */
};

/* 30 Nov 86 -- error in third coefficient fixed */

double zeta(double x, double q) {
  int i;
  double a, b, k, s, t, w;

  if (x == 1.0)
    goto retinf;

  if (x < 1.0) {
  domerr:
    // sf_error("zeta", SF_ERROR_DOMAIN, NULL);
    printf("ERROR");
    // return (NPY_NAN);
    return (NAN);
  }

  if (q <= 0.0) {
    if (q == floor(q)) {
      // sf_error("zeta", SF_ERROR_SINGULAR, NULL);
      printf("ERROR");
    retinf:
      // return (NPY_INFINITY);
      return (INFINITY);
    }
    if (x != floor(x))
      goto domerr; /* because q^-x not defined */
  }

  /* Asymptotic expansion
   * https://dlmf.nist.gov/25.11#E43
   */
  if (q > 1e8) {
    return (1 / (x - 1) + 1 / (2 * q)) * pow(q, 1 - x);
  }

  /* Euler-Maclaurin summation formula */

  /* Permit negative q but continue sum until n+q > +9 .
   * This case should be handled by a reflection formula.
   * If q<0 and x is an integer, there is a relation to
   * the polyGamma function.
   */
  s = pow(q, -x);
  a = q;
  i = 0;
  b = 0.0;
  while ((i < 9) || (a <= 9.0)) {
    i += 1;
    a += 1.0;
    b = pow(a, -x);
    s += b;
    if (fabs(b / s) < MACHEP)
      goto done;
  }

  w = a;
  s += b * w / (x - 1.0);
  s -= 0.5 * b;
  a = 1.0;
  k = 0.0;
  for (i = 0; i < 12; i++) {
    a *= x + k;
    b /= w;
    t = a * b / A_ZETA[i];
    s = s + t;
    t = fabs(t / s);
    if (t < MACHEP)
      goto done;
    k += 1.0;
    a *= x + k;
    b /= w;
    k += 1.0;
  }
done:
  return (s);
}
