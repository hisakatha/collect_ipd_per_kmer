#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
//#include <stdint.h>
#include <math.h>

// Collect IPD values by k-mer.
// k-mer is represented as a number with a radix of the size of chars,
// and the sum and the count of IPD values for a k-mer is return in an array at the k-mer index.
// Given a k-mer 5'-x_1 x_2 ... x_k-3',
// the k-mer index is T[x_1] * chars_size ^ (k-1) + T[x_2] * chars_size ^ (k-2) + ... + T[x_k],
// where, for example,
// x    | A C G T
// T[x] | 0 1 2 3
//
// coverage_threshold: IPD with coverage >= coverage_threshold will be used
// check_outside_coverage: whether to check coverage condition outside k-mer (1: true)
void collect_ipd_by_kmer(size_t const k, char const *chars, float const *tMeans, char **bases, size_t const dim,
        double *tMean_sum, double *tMean_sq_sum, double *tMean_log2_sum, double *tMean_log2_sq_sum,
        double *prediction_sum, double *prediction_sq_sum, double *prediction_log2_sum, double *prediction_log2_sq_sum, size_t *count, float const *modelPredictions,
        unsigned int const *coverage, unsigned int const coverage_threshold,
        size_t const outside_length, int const check_outside_coverage) {
    assert(dim % 2 == 0);
    assert(k < dim / 2);
    size_t chars_size = strlen(chars);
    size_t total_length = k + 2 * outside_length;
    // Prepare context holders for positive and negative strands
    int pos_context[k];
    int neg_context[k];
    int *context;
    // Indicate how many bases is required to collect k successive bases with a valid IPD.
    // Set to k if the current base is a null character, which means that no valid IPD is at the base
    int pos_state = k;
    int neg_state = k;
    int *state;
    // holders of tMean values for the surrounding region (useless?)
    //float pos_tMean_holder[k + 2 * outside_length];
    //float neg_tMean_holder[k + 2 * outside_length];
    for (size_t i = 0; i < dim; i++) {
        // Detect the current strand
        int isPositive = (i % 2 == 0);
        if(isPositive) {
            context = pos_context;
            state = &pos_state;
        } else {
            context = neg_context;
            state = &neg_state;
        }
        // Update context
        // Each context is oriented from 5' to 3' of the positive strand
        // For example, given a sequence
        // positive: 5'- ... x_1 x_2 ... x_k ... -3'
        // negative: 3'- ... y_1 y_2 ... y_k ... -5',
        // then PacBio HDF5 files contain data arrays (such as bases of this code) in the order of x_1 y_1 x_2 y_2 ..., and
        // pos_context: x_1 x_2 ... x_k
        // neg_context: y_1 y_2 ... y_k
        char *base_char = strchr(chars, bases[i][0]);
        unsigned int cur_coverage = coverage[i];
        if(base_char == NULL) {
            fprintf(stderr, "ERROR: Unexpected base was observed: %c\n", bases[i][0]);
            exit(EXIT_FAILURE);
        } else {
            for (size_t j = 0; j < k - 1; j++) {
                context[j] = context[j + 1];
            }
            context[k - 1] = base_char - chars;
            if (context[k - 1] == chars_size || cur_coverage < coverage_threshold) {
                //assert(bases[i][0] == '\0');
                *state = k;
            } else {
                *state = *state - 1;
            }
        }
        // Collect IPD using context
        if(*state <= 0) {
            // Reset to avoid negative overflow
            *state = 0;
            size_t sum_idx = 0;
            for (size_t j = 0; j < k; j++) {
                size_t context_idx = (isPositive) ? j : k - 1 - j;
                sum_idx = chars_size * sum_idx + context[context_idx];
            }
            sum_idx *= total_length;
            size_t tMean_idx_min = i - 2 * (k + outside_length - 1);
            size_t tMean_idx_max = i + 2 * outside_length;
            if(isPositive) {
                for (size_t tMean_idx = (tMean_idx_min < 0) ? 0 : tMean_idx_min; tMean_idx < dim && tMean_idx <= tMean_idx_max; tMean_idx += 2) {
                    double tMean = tMeans[tMean_idx];
                    double prediction = modelPredictions[tMean_idx];
                    if (tMean > 0.0 && (check_outside_coverage != 1 || coverage[tMean_idx] >= coverage_threshold)) {
                        tMean_sum[sum_idx] += tMean;
                        tMean_sq_sum[sum_idx] += tMean * tMean;
                        double tMean_log2 = log2(tMean);
                        tMean_log2_sum[sum_idx] += tMean_log2;
                        tMean_log2_sq_sum[sum_idx] += tMean_log2 * tMean_log2;
                        prediction_sum[sum_idx] += prediction;
                        prediction_sq_sum[sum_idx] += prediction * prediction;
                        double prediction_log2 = log2(prediction);
                        prediction_log2_sum[sum_idx] += prediction_log2;
                        prediction_log2_sq_sum[sum_idx] += prediction_log2 * prediction_log2;
                        count[sum_idx] += 1;
                    }
                    ++sum_idx;
                }
            } else {
                for (size_t tMean_idx = (tMean_idx_max >= dim) ? dim - 1 : tMean_idx_max; tMean_idx >= 0 && tMean_idx >= tMean_idx_min; tMean_idx -= 2) {
                    double tMean = tMeans[tMean_idx];
                    double prediction = modelPredictions[tMean_idx];
                    if (tMean > 0.0 && (check_outside_coverage != 1 || coverage[tMean_idx] >= coverage_threshold)) {
                        tMean_sum[sum_idx] += tMean;
                        tMean_sq_sum[sum_idx] += tMean * tMean;
                        double tMean_log2 = log2(tMean);
                        tMean_log2_sum[sum_idx] += tMean_log2;
                        tMean_log2_sq_sum[sum_idx] += tMean_log2 * tMean_log2;
                        prediction_sum[sum_idx] += prediction;
                        prediction_sq_sum[sum_idx] += prediction * prediction;
                        double prediction_log2 = log2(prediction);
                        prediction_log2_sum[sum_idx] += prediction_log2;
                        prediction_log2_sq_sum[sum_idx] += prediction_log2 * prediction_log2;
                        count[sum_idx] += 1;
                    }
                    ++sum_idx;
                }
            }
        }
    }
    return;
}
