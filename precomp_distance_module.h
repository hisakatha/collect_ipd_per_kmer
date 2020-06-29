#ifndef PRECOMP_DISTANCE_MODULE_H
#define PRECOMP_DISTANCE_MODULE_H

#ifdef __cplusplus
extern "C" {
#endif

    void collect_ipd_by_kmer(size_t k, size_t p, char const *chars, float const *tMeans, char **bases, size_t dim, double *sum, size_t *count, char const *dataset_output_path, float const *modelPredictions, char const *pacbio_model_output_path);
    size_t insert_base_into_context(size_t context, size_t chars_size, size_t k, size_t position, size_t base);
    double compute_base_distance(size_t position, size_t base1, size_t base2, size_t chars_size, size_t k, double const *ipd_sum, size_t const * count, double count_threshold_rate, size_t *min_valid_context_count);
    void fill_base_distance_matrix(size_t chars_size, size_t k, double const * ipd_sum, size_t const * count, double count_threshold_rate, double position_relative_weight, double *matrix);

#ifdef __cplusplus
}
#endif

#endif
