#ifndef POLEVL_HPP
#define POLEVL_HPP

double polevl(double x, double const* coef, int N_pol);
double p1evl(double x, const double coef[], int N_pol);
double ratevl(double x, const double num[], int M,
                            const double denom[], int N_pol);

#endif