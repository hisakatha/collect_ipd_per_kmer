#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

// Collect IPD values by k-mer.
// k-mer is represented as a number with a radix of the size of chars,
// and the sum and the count of IPD values for a k-mer is return in an array at the k-mer index.
// Given a k-mer 5'-x_1 x_2 ... x_k-3',
// the k-mer index is T[x_1] * chars_size ^ (k-1) + T[x_2] * chars_size ^ (k-2) + ... + T[x_k],
// where, for example,
// x    | A C G T
// T[x] | 0 1 2 3
void collect_ipd_by_kmer(size_t k, size_t p, char const *chars, float const *tMeans, char **bases, size_t dim, double *sum, size_t *count, char const *dataset_output_path, float const *modelPredictions, char const *pacbio_model_output_path) {
    assert(dim % 2 == 0);
    assert(k < dim / 2);
    assert(p <= k);
    size_t chars_size = strlen(chars);
    FILE *dataset_out = NULL;
    if(dataset_output_path != NULL){
        dataset_out = fopen(dataset_output_path, "a");
        if(dataset_out == NULL){
            fprintf(stderr, "ERROR: Cannot open file: %s\n", dataset_output_path); exit(EXIT_FAILURE);
        }
        // Set a buffer size for speedup
        if(setvbuf(dataset_out, NULL, _IOFBF, 16 * 1024 * 1024) != 0){
            fprintf(stderr, "WARNING: Setting a buffer size for data set failed\n");
            // Reset the buffer size to BUFSIZ defined in stdio.h
            if(setvbuf(dataset_out, NULL, _IOFBF, BUFSIZ) != 0){
                fprintf(stderr, "ERROR: Cannot set a buffer size for data set\n"); exit(EXIT_FAILURE);
            }
        }
    }
    FILE *pacbio_model_out = NULL;
    if(pacbio_model_output_path != NULL){
        pacbio_model_out = fopen(pacbio_model_output_path, "a");
        if(pacbio_model_out == NULL){
            fprintf(stderr, "ERROR: Cannot open file: %s\n", pacbio_model_output_path); exit(EXIT_FAILURE);
        }
        // Set a buffer size for speedup
        if(setvbuf(pacbio_model_out, NULL, _IOFBF, 16 * 1024 * 1024) != 0){
            fprintf(stderr, "WARNING: Setting a buffer size for pacbio model output failed\n");
            // Reset the buffer size to BUFSIZ defined in stdio.h
            if(setvbuf(pacbio_model_out, NULL, _IOFBF, BUFSIZ) != 0){
                fprintf(stderr, "ERROR: Cannot set a buffer size for pacbio model output\n"); exit(EXIT_FAILURE);
            }
        }
    }
    // Prepare context holders for positive and negative strands
    int pos_context[k];
    int neg_context[k];
    int *context;
    // Indicate how many bases is required to collect k successive bases with a valid IPD.
    // Set to k if the current base is a null character, which means that no valid IPD is at the base
    int pos_state = k;
    int neg_state = k;
    int *state;
    for (size_t i = 0; i < dim; i++) {
        // Detect the current strand
        if(i % 2 == 0) {
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
        if(base_char == NULL) {
            fprintf(stderr, "ERROR: Unexpected base was observed: %c\n", bases[i][0]);
            exit(EXIT_FAILURE);
        } else {
            for (size_t j = 0; j < k - 1; j++) {
                context[j] = context[j + 1];
            }
            context[k - 1] = base_char - chars;
            if (context[k - 1] == chars_size) {
                assert(bases[i][0] == '\0');
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
                size_t context_idx = (i % 2 == 0) ? j : k - 1 - j;
                sum_idx = chars_size * sum_idx + context[context_idx];
            }
            // Get the IPD of the p-th base in the k-mer whose last base is at 0-based position i
            // If it is on the positive strand, move back (k-1) - (p-1); if not, move back p-1
            size_t tMean_idx = i - 2 * (i % 2 == 0 ? k - p : p - 1);
            float tMean = tMeans[tMean_idx];
            sum[sum_idx] += tMean;
            count[sum_idx] += 1;
            if(dataset_out != NULL){
                fprintf(dataset_out, "%.9g 1:%zu\n", tMean, sum_idx);
            }
            if(pacbio_model_out != NULL){
                fprintf(pacbio_model_out, "%.9g\n", modelPredictions[tMean_idx]);
            }
        }
    }
    if(dataset_out != NULL){ fclose(dataset_out); }
    if(pacbio_model_out != NULL){ fclose(pacbio_model_out); }
    return;
}

// Return a context represented as a number with a radix of chars_size
// The resulting context is a context of input context into which a base is inserted
// position: 0-based index which indicates where to insert a base
size_t insert_base_into_context(size_t context, size_t chars_size, size_t k, size_t position, size_t base){
    assert(position < k);
    size_t insert_coef = (size_t)(pow(chars_size, k - 1 - position) + 0.5);
    return ((context / insert_coef) * insert_coef * chars_size) + base * insert_coef + context % insert_coef;
}

// If there are enough k-mers (according to a threshold), then return a distance between bases;
// otherwise, return -1.0
double compute_base_distance(size_t position, size_t base1, size_t base2, size_t chars_size, size_t k, double const * ipd_sum, size_t const * count, double count_threshold_rate, size_t *min_valid_context_count){
    // Scan (k-1)-mer space
    size_t max_context = (size_t)(pow(chars_size, k - 1) + 0.5);
    size_t valid_context_count = 0;
    double current_sum = 0.0;
    for(size_t orig_context = 0; orig_context < max_context; orig_context++){
        size_t context1 = insert_base_into_context(orig_context, chars_size, k, position, base1);
        size_t context2 = insert_base_into_context(orig_context, chars_size, k, position, base2);
        size_t count1 = count[context1];
        size_t count2 = count[context2];
        // TODO: This threshold should be parameterized?
        if(count1 > 0 && count2 > 0){
            valid_context_count++;
            double d = ipd_sum[context1] / count1 - ipd_sum[context2] / count2;
            current_sum += d * d;
        }
    }
    double distance = -1.0;
    if(valid_context_count > count_threshold_rate * (max_context - 1)){
        distance = sqrt(current_sum / valid_context_count);
    }
    if(*min_valid_context_count > valid_context_count){
        *min_valid_context_count = valid_context_count;
    }
    return distance;
}

// Fill a distance matrix between bases,
// and complement a missing distance using other elements if possible
void fill_base_distance_matrix(size_t chars_size, size_t k, double const * ipd_sum, size_t const * count, double count_threshold_rate, double position_relative_weight, double *matrix){
    // Fill a distance matrix without complementing
    size_t min_valid_context_count = SIZE_MAX;
    for(size_t position = 0; position < k; position++){
        for(size_t base1 = 0; base1 < chars_size; base1++){
            for(size_t base2 = base1; base2 < chars_size; base2++){
                size_t idx = position * chars_size * chars_size + base1 * chars_size + base2;
                if(base1 == base2){
                    matrix[idx] = 0.0;
                } else {
                    matrix[idx] = compute_base_distance(position, base1, base2, chars_size, k, ipd_sum, count, count_threshold_rate, &min_valid_context_count);
                    size_t counter_idx = position * chars_size * chars_size + base2 * chars_size + base1;
                    matrix[counter_idx] = matrix[idx];
                }
            }
        }
    }
    fprintf(stderr, "INFO: min_valid_context_count = %zu, min_count_rate = %g\n", min_valid_context_count, min_valid_context_count / (pow(chars_size, k - 1) - 1.0));
    // Complement missing distances
    size_t num_count_threshold_failure = 0;
    for(size_t position = 0; position < k; position++){
        for(size_t base1 = 0; base1 < chars_size; base1++){
            for(size_t base2 = base1 + 1; base2 < chars_size; base2++){
                size_t idx = position * chars_size * chars_size + base1 * chars_size + base2;
                if(matrix[idx] < 0.0){
                    num_count_threshold_failure++;
                    // Complement using position-fixed contexts
                    size_t count_fixed_position = 0;
                    double sum_fixed_position = 0.0;
                    for(size_t base1another = 0; base1another < chars_size; base1another++){
                        for(size_t base2another = base1another + 1; base2another < chars_size; base2another++){
                            size_t idx_fixed_position = position * chars_size * chars_size + base1another * chars_size + base2another;
                            if(matrix[idx_fixed_position] > 0.0){
                                count_fixed_position++;
                                sum_fixed_position += matrix[idx_fixed_position];
                            }
                        }
                    }
                    // Complement using bases-fixed contexts
                    size_t count_fixed_bases = 0;
                    double sum_fixed_bases = 0.0;
                    for(size_t position_another = 0; position_another < k; position_another++){
                        size_t idx_fixed_bases = position_another * chars_size * chars_size + base1 * chars_size + base2;
                        if(matrix[idx_fixed_bases] > 0.0){
                            count_fixed_bases++;
                            sum_fixed_bases += matrix[idx_fixed_bases];
                        }
                    }
                    // Combine above complements
                    if(count_fixed_position > 0 && count_fixed_bases > 0){
                        matrix[idx] = position_relative_weight * sum_fixed_position / count_fixed_position + (1 - position_relative_weight) * sum_fixed_bases / count_fixed_bases;
                    } else {
                        fprintf(stderr, "ERROR: No support for complementing a distance of position: %zu, base1: %zu, base2: %zu\n", position, base1, base2);
                        matrix[idx] = - INFINITY;
                    }
                    size_t counter_idx = position * chars_size * chars_size + base2 * chars_size + base1;
                    matrix[counter_idx] = matrix[idx];
                }
            }
        }
    }
    fprintf(stderr, "INFO: # element in the upper triangle of a base distance to be complemented: %zu\n", num_count_threshold_failure);
    return;
}
