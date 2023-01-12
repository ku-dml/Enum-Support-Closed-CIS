#ifndef PVALUE_HPP
#define PVALUE_HPP

#include "../clost/define.h"

class _Stat {
    public:
        // Initialize class
        _Stat(Tool T, Graph G);
        // Compute isf(alpha / k) where isf is the inverse of survive function
        double inverse_threshold(double alpha, int k);
        // Compute p-value
        double survival_function(double x);
        // Compute isf(p-value of a given S)
        double p_value(OwnStack S);
        // Compute isf(minimal p-value of a given S)
        double minimal_p_value(OwnStack S);
        // Compute isf(envelope of a given S)
        double envelope(OwnStack S);
    private:
        // only use in order to compute envelope
        double minimal_p_value_inner(const std::vector<double>& xs);
        // given data
        int J;
        Tool tool;
        Graph graph;
        // use to avoid unnecessary heap allocations
        itemset I_buffer;
        std::vector<double> xs_buffer;
        std::vector<std::tuple<double, int>> betaLs_buffer;
        std::vector<std::tuple<double, int>> betaRs_buffer;
        std::vector<double> xStarsL_buffer;
        std::vector<double> xStarsR_buffer;
};
typedef _Stat *Stat;

#endif