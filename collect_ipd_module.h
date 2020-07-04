#ifndef COLLECT_IPD_MODULE_H
#define COLLECT_IPD_MODULE_H

#ifdef __cplusplus
extern "C" {
#endif

    void collect_ipd_by_kmer(size_t const k, char const *chars, float const *tMeans, char **bases, size_t const dim,
        double *tMean_sum, double *tMean_sq_sum, double *tMean_log2_sum, double *tMean_log2_sq_sum,
        double *prediction_sum, double *prediction_sq_sum, double *prediction_log2_sum, double *prediction_log2_sq_sum, size_t *count, float const *modelPredictions,
        unsigned int const *coverage, unsigned int const coverage_threshold,
        size_t const outside_length, int const check_outside_coverage);

#ifdef __cplusplus
}
#endif

#endif
