#include <stdint.h>
#include <cmath>
#include <CppUTest/CommandLineTestRunner.h>
#include "collect_ipd_module.h"

TEST_GROUP(kmer_ipd)
{
    // Tests for
    // void collect_ipd_by_kmer(size_t const k, char const *chars, float const *tMeans, char **bases, size_t const dim,
    //    double *tMean_sum, double *tMean_sq_sum, double *tMean_log2_sum, double *tMean_log2_sq_sum,
    //    double *prediction_sum, double *prediction_sq_sum, double *prediction_log2_sum, double *prediction_log2_sq_sum, size_t *count, float const *modelPredictions,
    //    unsigned int const *coverage, unsigned int const coverage_threshold,
    //    size_t const outside_length, int const check_outside_coverage);
    char *chars = (char *)"ACGT";
    size_t dim = 10;
    size_t coverage_threshold = 25;
    double tolerance = 0.0001;

    // Input arrays correspond to {pos1, neg1, pos2, neg2, ..., pos_n, neg_n}, that is,
    // positive: AACGC, negative: GC0TT
    const char *bases[10]      = {"A", "T", "A", "T", "C", "\0", "G", "C", "C", "G"};
    float tMeans[10]           = {1.5, 2.2, 3.1, 4.9, 5.5, 6.3, 7.2, 8.2, 9.1, 10.3};
    float modelPredictions[10] = {1.4, 2.1, 3.0, 4.8, 5.4, 6.2, 7.1, 8.1, 9.0, 10.2};
    unsigned int coverage[10]  = { 30,  20,  30,  20,  30,  20,  30,  40,  40,   40};
};

TEST(kmer_ipd, k1l0)
{
    size_t k = 1;
    size_t outside_length = 0;
    size_t total_length = k + 2 * outside_length;
    int check_outside_coverage = 0;
    size_t kmers_size = (size_t)(std::pow(4, k) + 0.5);
    CHECK_EQUAL(kmers_size, 4);
    size_t array_size = kmers_size * total_length;
    double tMean_sum[array_size] = {0.0};
    double tMean_sq_sum[array_size] = {0.0};
    double tMean_log2_sum[array_size] = {0.0};
    double tMean_log2_sq_sum[array_size] = {0.0};
    double prediction_sum[array_size] = {0.0};
    double prediction_sq_sum[array_size] = {0.0};
    double prediction_log2_sum[array_size] = {0.0};
    double prediction_log2_sq_sum[array_size] = {0.0};
    size_t count[array_size] = {0};
    collect_ipd_by_kmer(k, chars, tMeans, (char **)bases, dim, tMean_sum, tMean_sq_sum, tMean_log2_sum, tMean_log2_sq_sum,
            prediction_sum, prediction_sq_sum, prediction_log2_sum, prediction_log2_sq_sum, count, modelPredictions, coverage, coverage_threshold, outside_length, check_outside_coverage);
    CHECK_EQUAL(1.5f + 3.1f, (float)tMean_sum[0]);
    CHECK_EQUAL(5.5f + 8.2f + 9.1f, (float)tMean_sum[1]);
    CHECK_EQUAL(7.2f + 10.3f, (float)tMean_sum[2]);
    CHECK_EQUAL(0.0, (float)tMean_sum[3]);
    printf("all_double: %.17g, float_to_double: %.17g, all_float, %.17g, actual: %.17g\n", 1.5*1.5+3.1*3.1, (double)1.5f*(double)1.5f + (double)3.1f*(double)3.1f, 1.5f*1.5f+3.1f*3.1f, tMean_sq_sum[0]);
    DOUBLES_EQUAL(1.5*1.5 + 3.1*3.1, tMean_sq_sum[0], tolerance);
    DOUBLES_EQUAL(5.5*5.5 + 8.2*8.2 + 9.1*9.1, tMean_sq_sum[1], tolerance);
    DOUBLES_EQUAL(7.2*7.2 + 10.3*10.3, tMean_sq_sum[2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[3], tolerance);
    DOUBLES_EQUAL(log2(1.5) + log2(3.1), tMean_log2_sum[0], tolerance);
    DOUBLES_EQUAL(log2(5.5) + log2(8.2) + log2(9.1), tMean_log2_sum[1], tolerance);
    DOUBLES_EQUAL(log2(7.2) + log2(10.3), tMean_log2_sum[2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_log2_sum[3], tolerance);
    DOUBLES_EQUAL(log2(1.5)*log2(1.5) + log2(3.1)*log2(3.1), tMean_log2_sq_sum[0], tolerance);
    DOUBLES_EQUAL(log2(5.5)*log2(5.5) + log2(8.2)*log2(8.2) + log2(9.1)*log2(9.1), tMean_log2_sq_sum[1], tolerance);
    DOUBLES_EQUAL(log2(7.2)*log2(7.2) + log2(10.3)*log2(10.3), tMean_log2_sq_sum[2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_log2_sum[3], tolerance);
    DOUBLES_EQUAL(1.4 + 3.0, prediction_sum[0], tolerance);
    DOUBLES_EQUAL(5.4 + 8.1 + 9.0, prediction_sum[1], tolerance);
    DOUBLES_EQUAL(7.1 + 10.2, prediction_sum[2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[3], tolerance);
    DOUBLES_EQUAL(1.4*1.4 + 3.0*3.0, prediction_sq_sum[0], tolerance);
    DOUBLES_EQUAL(5.4*5.4 + 8.1*8.1 + 9.0*9.0, prediction_sq_sum[1], tolerance);
    DOUBLES_EQUAL(7.1*7.1 + 10.2*10.2, prediction_sq_sum[2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[3], tolerance);
    DOUBLES_EQUAL(log2(1.4) + log2(3.0), prediction_log2_sum[0], tolerance);
    DOUBLES_EQUAL(log2(5.4) + log2(8.1) + log2(9.0), prediction_log2_sum[1], tolerance);
    DOUBLES_EQUAL(log2(7.1) + log2(10.2), prediction_log2_sum[2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_log2_sum[3], tolerance);
    DOUBLES_EQUAL(log2(1.4)*log2(1.4) + log2(3.0)*log2(3.0), prediction_log2_sq_sum[0], tolerance);
    DOUBLES_EQUAL(log2(5.4)*log2(5.4) + log2(8.1)*log2(8.1) + log2(9.0)*log2(9.0), prediction_log2_sq_sum[1], tolerance);
    DOUBLES_EQUAL(log2(7.1)*log2(7.1) + log2(10.2)*log2(10.2), prediction_log2_sq_sum[2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_log2_sum[3], tolerance);
    CHECK_EQUAL(2, count[0]);
    CHECK_EQUAL(3, count[1]);
    CHECK_EQUAL(2, count[2]);
    CHECK_EQUAL(0, count[3]);
}

TEST(kmer_ipd, k2l1)
{
    size_t k = 2;
    size_t outside_length = 1;
    size_t total_length = k + 2 * outside_length;
    int check_outside_coverage = 1;
    size_t kmers_size = (size_t)(std::pow(4, k) + 0.5);
    CHECK_EQUAL(16, kmers_size);
    size_t array_size = kmers_size * total_length;
    double tMean_sum[array_size] = {0.0};
    double tMean_sq_sum[array_size] = {0.0};
    double tMean_log2_sum[array_size] = {0.0};
    double tMean_log2_sq_sum[array_size] = {0.0};
    double prediction_sum[array_size] = {0.0};
    double prediction_sq_sum[array_size] = {0.0};
    double prediction_log2_sum[array_size] = {0.0};
    double prediction_log2_sq_sum[array_size] = {0.0};
    size_t count[array_size] = {0};
    collect_ipd_by_kmer(k, chars, tMeans, (char **)bases, dim, tMean_sum, tMean_sq_sum, tMean_log2_sum, tMean_log2_sq_sum,
            prediction_sum, prediction_sq_sum, prediction_log2_sum, prediction_log2_sq_sum, count, modelPredictions, coverage, coverage_threshold, outside_length, check_outside_coverage);
    // AA
    DOUBLES_EQUAL(0.0, tMean_sum[0 * total_length + 0], tolerance);
    DOUBLES_EQUAL(1.5, tMean_sum[0 * total_length + 1], tolerance);
    DOUBLES_EQUAL(3.1, tMean_sum[0 * total_length + 2], tolerance);
    DOUBLES_EQUAL(5.5, tMean_sum[0 * total_length + 3], tolerance);
    // AC
    DOUBLES_EQUAL(1.5, tMean_sum[1 * total_length + 0], tolerance);
    DOUBLES_EQUAL(3.1, tMean_sum[1 * total_length + 1], tolerance);
    DOUBLES_EQUAL(5.5, tMean_sum[1 * total_length + 2], tolerance);
    DOUBLES_EQUAL(7.2, tMean_sum[1 * total_length + 3], tolerance);
    // AG
    DOUBLES_EQUAL(0.0, tMean_sum[2 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[2 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[2 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[2 * total_length + 3], tolerance);
    // AT
    DOUBLES_EQUAL(0.0, tMean_sum[3 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[3 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[3 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[3 * total_length + 3], tolerance);
    // CA
    DOUBLES_EQUAL(0.0, tMean_sum[4 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[4 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[4 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[4 * total_length + 3], tolerance);
    // CC
    DOUBLES_EQUAL(0.0, tMean_sum[5 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[5 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[5 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[5 * total_length + 3], tolerance);
    // CG
    DOUBLES_EQUAL(3.1, tMean_sum[6 * total_length + 0], tolerance);
    DOUBLES_EQUAL(5.5, tMean_sum[6 * total_length + 1], tolerance);
    DOUBLES_EQUAL(7.2, tMean_sum[6 * total_length + 2], tolerance);
    DOUBLES_EQUAL(9.1, tMean_sum[6 * total_length + 3], tolerance);
    // CT
    DOUBLES_EQUAL(0.0, tMean_sum[7 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[7 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[7 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[7 * total_length + 3], tolerance);
    // GA
    DOUBLES_EQUAL(0.0, tMean_sum[8 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[8 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[8 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[8 * total_length + 3], tolerance);
    // GC
    DOUBLES_EQUAL(5.5, tMean_sum[9 * total_length + 0], tolerance);
    DOUBLES_EQUAL(7.2 + 10.3, tMean_sum[9 * total_length + 1], tolerance);
    DOUBLES_EQUAL(9.1 + 8.2, tMean_sum[9 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[9 * total_length + 3], tolerance);
    // GG
    DOUBLES_EQUAL(0.0, tMean_sum[10 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[10 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[10 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[10 * total_length + 3], tolerance);
    // GT
    DOUBLES_EQUAL(0.0, tMean_sum[11 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[11 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[11 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[11 * total_length + 3], tolerance);
    // TA
    DOUBLES_EQUAL(0.0, tMean_sum[12 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[12 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[12 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[12 * total_length + 3], tolerance);
    // TC
    DOUBLES_EQUAL(0.0, tMean_sum[13 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[13 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[13 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[13 * total_length + 3], tolerance);
    // TG
    DOUBLES_EQUAL(0.0, tMean_sum[14 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[14 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[14 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[14 * total_length + 3], tolerance);
    // TT
    DOUBLES_EQUAL(0.0, tMean_sum[15 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[15 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[15 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sum[15 * total_length + 3], tolerance);

    // AA
    DOUBLES_EQUAL(0.0, tMean_sq_sum[0 * total_length + 0], tolerance);
    DOUBLES_EQUAL(1.5 * 1.5, tMean_sq_sum[0 * total_length + 1], tolerance);
    DOUBLES_EQUAL(3.1 * 3.1, tMean_sq_sum[0 * total_length + 2], tolerance);
    DOUBLES_EQUAL(5.5 * 5.5, tMean_sq_sum[0 * total_length + 3], tolerance);
    // AC
    DOUBLES_EQUAL(1.5 * 1.5, tMean_sq_sum[1 * total_length + 0], tolerance);
    DOUBLES_EQUAL(3.1 * 3.1, tMean_sq_sum[1 * total_length + 1], tolerance);
    DOUBLES_EQUAL(5.5 * 5.5, tMean_sq_sum[1 * total_length + 2], tolerance);
    DOUBLES_EQUAL(7.2 * 7.2, tMean_sq_sum[1 * total_length + 3], tolerance);
    // AG
    DOUBLES_EQUAL(0.0, tMean_sq_sum[2 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[2 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[2 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[2 * total_length + 3], tolerance);
    // AT
    DOUBLES_EQUAL(0.0, tMean_sq_sum[3 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[3 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[3 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[3 * total_length + 3], tolerance);
    // CA
    DOUBLES_EQUAL(0.0, tMean_sq_sum[4 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[4 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[4 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[4 * total_length + 3], tolerance);
    // CC
    DOUBLES_EQUAL(0.0, tMean_sq_sum[5 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[5 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[5 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[5 * total_length + 3], tolerance);
    // CG
    DOUBLES_EQUAL(3.1 * 3.1, tMean_sq_sum[6 * total_length + 0], tolerance);
    DOUBLES_EQUAL(5.5 * 5.5, tMean_sq_sum[6 * total_length + 1], tolerance);
    DOUBLES_EQUAL(7.2 * 7.2, tMean_sq_sum[6 * total_length + 2], tolerance);
    DOUBLES_EQUAL(9.1 * 9.1, tMean_sq_sum[6 * total_length + 3], tolerance);
    // CT
    DOUBLES_EQUAL(0.0, tMean_sq_sum[7 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[7 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[7 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[7 * total_length + 3], tolerance);
    // GA
    DOUBLES_EQUAL(0.0, tMean_sq_sum[8 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[8 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[8 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[8 * total_length + 3], tolerance);
    // GC
    DOUBLES_EQUAL(5.5 * 5.5, tMean_sq_sum[9 * total_length + 0], tolerance);
    DOUBLES_EQUAL(7.2 * 7.2 + 10.3 * 10.3, tMean_sq_sum[9 * total_length + 1], tolerance);
    DOUBLES_EQUAL(9.1 * 9.1 + 8.2 * 8.2, tMean_sq_sum[9 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[9 * total_length + 3], tolerance);
    // GG
    DOUBLES_EQUAL(0.0, tMean_sq_sum[10 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[10 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[10 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[10 * total_length + 3], tolerance);
    // GT
    DOUBLES_EQUAL(0.0, tMean_sq_sum[11 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[11 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[11 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[11 * total_length + 3], tolerance);
    // TA
    DOUBLES_EQUAL(0.0, tMean_sq_sum[12 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[12 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[12 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[12 * total_length + 3], tolerance);
    // TC
    DOUBLES_EQUAL(0.0, tMean_sq_sum[13 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[13 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[13 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[13 * total_length + 3], tolerance);
    // TG
    DOUBLES_EQUAL(0.0, tMean_sq_sum[14 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[14 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[14 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[14 * total_length + 3], tolerance);
    // TT
    DOUBLES_EQUAL(0.0, tMean_sq_sum[15 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[15 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[15 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, tMean_sq_sum[15 * total_length + 3], tolerance);

    // AA
    DOUBLES_EQUAL(0.0, prediction_sum[0 * total_length + 0], tolerance);
    DOUBLES_EQUAL(1.4, prediction_sum[0 * total_length + 1], tolerance);
    DOUBLES_EQUAL(3.0, prediction_sum[0 * total_length + 2], tolerance);
    DOUBLES_EQUAL(5.4, prediction_sum[0 * total_length + 3], tolerance);
    // AC
    DOUBLES_EQUAL(1.4, prediction_sum[1 * total_length + 0], tolerance);
    DOUBLES_EQUAL(3.0, prediction_sum[1 * total_length + 1], tolerance);
    DOUBLES_EQUAL(5.4, prediction_sum[1 * total_length + 2], tolerance);
    DOUBLES_EQUAL(7.1, prediction_sum[1 * total_length + 3], tolerance);
    // AG
    DOUBLES_EQUAL(0.0, prediction_sum[2 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[2 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[2 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[2 * total_length + 3], tolerance);
    // AT
    DOUBLES_EQUAL(0.0, prediction_sum[3 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[3 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[3 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[3 * total_length + 3], tolerance);
    // CA
    DOUBLES_EQUAL(0.0, prediction_sum[4 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[4 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[4 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[4 * total_length + 3], tolerance);
    // CC
    DOUBLES_EQUAL(0.0, prediction_sum[5 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[5 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[5 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[5 * total_length + 3], tolerance);
    // CG
    DOUBLES_EQUAL(3.0, prediction_sum[6 * total_length + 0], tolerance);
    DOUBLES_EQUAL(5.4, prediction_sum[6 * total_length + 1], tolerance);
    DOUBLES_EQUAL(7.1, prediction_sum[6 * total_length + 2], tolerance);
    DOUBLES_EQUAL(9.0, prediction_sum[6 * total_length + 3], tolerance);
    // CT
    DOUBLES_EQUAL(0.0, prediction_sum[7 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[7 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[7 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[7 * total_length + 3], tolerance);
    // GA
    DOUBLES_EQUAL(0.0, prediction_sum[8 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[8 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[8 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[8 * total_length + 3], tolerance);
    // GC
    DOUBLES_EQUAL(5.4, prediction_sum[9 * total_length + 0], tolerance);
    DOUBLES_EQUAL(7.1 + 10.2, prediction_sum[9 * total_length + 1], tolerance);
    DOUBLES_EQUAL(9.0 + 8.1, prediction_sum[9 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[9 * total_length + 3], tolerance);
    // GG
    DOUBLES_EQUAL(0.0, prediction_sum[10 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[10 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[10 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[10 * total_length + 3], tolerance);
    // GT
    DOUBLES_EQUAL(0.0, prediction_sum[11 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[11 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[11 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[11 * total_length + 3], tolerance);
    // TA
    DOUBLES_EQUAL(0.0, prediction_sum[12 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[12 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[12 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[12 * total_length + 3], tolerance);
    // TC
    DOUBLES_EQUAL(0.0, prediction_sum[13 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[13 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[13 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[13 * total_length + 3], tolerance);
    // TG
    DOUBLES_EQUAL(0.0, prediction_sum[14 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[14 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[14 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[14 * total_length + 3], tolerance);
    // TT
    DOUBLES_EQUAL(0.0, prediction_sum[15 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[15 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[15 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sum[15 * total_length + 3], tolerance);

    // AA
    DOUBLES_EQUAL(0.0, prediction_sq_sum[0 * total_length + 0], tolerance);
    DOUBLES_EQUAL(1.4 * 1.4, prediction_sq_sum[0 * total_length + 1], tolerance);
    DOUBLES_EQUAL(3.0 * 3.0, prediction_sq_sum[0 * total_length + 2], tolerance);
    DOUBLES_EQUAL(5.4 * 5.4, prediction_sq_sum[0 * total_length + 3], tolerance);
    // AC
    DOUBLES_EQUAL(1.4 * 1.4, prediction_sq_sum[1 * total_length + 0], tolerance);
    DOUBLES_EQUAL(3.0 * 3.0, prediction_sq_sum[1 * total_length + 1], tolerance);
    DOUBLES_EQUAL(5.4 * 5.4, prediction_sq_sum[1 * total_length + 2], tolerance);
    DOUBLES_EQUAL(7.1 * 7.1, prediction_sq_sum[1 * total_length + 3], tolerance);
    // AG
    DOUBLES_EQUAL(0.0, prediction_sq_sum[2 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[2 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[2 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[2 * total_length + 3], tolerance);
    // AT
    DOUBLES_EQUAL(0.0, prediction_sq_sum[3 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[3 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[3 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[3 * total_length + 3], tolerance);
    // CA
    DOUBLES_EQUAL(0.0, prediction_sq_sum[4 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[4 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[4 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[4 * total_length + 3], tolerance);
    // CC
    DOUBLES_EQUAL(0.0, prediction_sq_sum[5 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[5 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[5 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[5 * total_length + 3], tolerance);
    // CG
    DOUBLES_EQUAL(3.0 * 3.0, prediction_sq_sum[6 * total_length + 0], tolerance);
    DOUBLES_EQUAL(5.4 * 5.4, prediction_sq_sum[6 * total_length + 1], tolerance);
    DOUBLES_EQUAL(7.1 * 7.1, prediction_sq_sum[6 * total_length + 2], tolerance);
    DOUBLES_EQUAL(9.0 * 9.0, prediction_sq_sum[6 * total_length + 3], tolerance);
    // CT
    DOUBLES_EQUAL(0.0, prediction_sq_sum[7 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[7 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[7 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[7 * total_length + 3], tolerance);
    // GA
    DOUBLES_EQUAL(0.0, prediction_sq_sum[8 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[8 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[8 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[8 * total_length + 3], tolerance);
    // GC
    DOUBLES_EQUAL(5.4 * 5.4, prediction_sq_sum[9 * total_length + 0], tolerance);
    DOUBLES_EQUAL(7.1 * 7.1 + 10.2 * 10.2, prediction_sq_sum[9 * total_length + 1], tolerance);
    DOUBLES_EQUAL(9.0 * 9.0 + 8.1 * 8.1, prediction_sq_sum[9 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[9 * total_length + 3], tolerance);
    // GG
    DOUBLES_EQUAL(0.0, prediction_sq_sum[10 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[10 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[10 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[10 * total_length + 3], tolerance);
    // GT
    DOUBLES_EQUAL(0.0, prediction_sq_sum[11 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[11 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[11 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[11 * total_length + 3], tolerance);
    // TA
    DOUBLES_EQUAL(0.0, prediction_sq_sum[12 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[12 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[12 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[12 * total_length + 3], tolerance);
    // TC
    DOUBLES_EQUAL(0.0, prediction_sq_sum[13 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[13 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[13 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[13 * total_length + 3], tolerance);
    // TG
    DOUBLES_EQUAL(0.0, prediction_sq_sum[14 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[14 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[14 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[14 * total_length + 3], tolerance);
    // TT
    DOUBLES_EQUAL(0.0, prediction_sq_sum[15 * total_length + 0], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[15 * total_length + 1], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[15 * total_length + 2], tolerance);
    DOUBLES_EQUAL(0.0, prediction_sq_sum[15 * total_length + 3], tolerance);

    // AA
    CHECK_EQUAL(0, count[0 * total_length + 0]);
    CHECK_EQUAL(1, count[0 * total_length + 1]);
    CHECK_EQUAL(1, count[0 * total_length + 2]);
    CHECK_EQUAL(1, count[0 * total_length + 3]);
    // AC
    CHECK_EQUAL(1, count[1 * total_length + 0]);
    CHECK_EQUAL(1, count[1 * total_length + 1]);
    CHECK_EQUAL(1, count[1 * total_length + 2]);
    CHECK_EQUAL(1, count[1 * total_length + 3]);
    // AG
    CHECK_EQUAL(0, count[2 * total_length + 0]);
    CHECK_EQUAL(0, count[2 * total_length + 1]);
    CHECK_EQUAL(0, count[2 * total_length + 2]);
    CHECK_EQUAL(0, count[2 * total_length + 3]);
    // AT
    CHECK_EQUAL(0, count[3 * total_length + 0]);
    CHECK_EQUAL(0, count[3 * total_length + 1]);
    CHECK_EQUAL(0, count[3 * total_length + 2]);
    CHECK_EQUAL(0, count[3 * total_length + 3]);
    // CA
    CHECK_EQUAL(0, count[4 * total_length + 0]);
    CHECK_EQUAL(0, count[4 * total_length + 1]);
    CHECK_EQUAL(0, count[4 * total_length + 2]);
    CHECK_EQUAL(0, count[4 * total_length + 3]);
    // CC
    CHECK_EQUAL(0, count[5 * total_length + 0]);
    CHECK_EQUAL(0, count[5 * total_length + 1]);
    CHECK_EQUAL(0, count[5 * total_length + 2]);
    CHECK_EQUAL(0, count[5 * total_length + 3]);
    // CG
    CHECK_EQUAL(1, count[6 * total_length + 0]);
    CHECK_EQUAL(1, count[6 * total_length + 1]);
    CHECK_EQUAL(1, count[6 * total_length + 2]);
    CHECK_EQUAL(1, count[6 * total_length + 3]);
    // CT
    CHECK_EQUAL(0, count[7 * total_length + 0]);
    CHECK_EQUAL(0, count[7 * total_length + 1]);
    CHECK_EQUAL(0, count[7 * total_length + 2]);
    CHECK_EQUAL(0, count[7 * total_length + 3]);
    // GA
    CHECK_EQUAL(0, count[8 * total_length + 0]);
    CHECK_EQUAL(0, count[8 * total_length + 1]);
    CHECK_EQUAL(0, count[8 * total_length + 2]);
    CHECK_EQUAL(0, count[8 * total_length + 3]);
    // GC
    CHECK_EQUAL(1, count[9 * total_length + 0]);
    CHECK_EQUAL(2, count[9 * total_length + 1]);
    CHECK_EQUAL(2, count[9 * total_length + 2]);
    CHECK_EQUAL(0, count[9 * total_length + 3]);
    // GG
    CHECK_EQUAL(0, count[10 * total_length + 0]);
    CHECK_EQUAL(0, count[10 * total_length + 1]);
    CHECK_EQUAL(0, count[10 * total_length + 2]);
    CHECK_EQUAL(0, count[10 * total_length + 3]);
    // GT
    CHECK_EQUAL(0, count[11 * total_length + 0]);
    CHECK_EQUAL(0, count[11 * total_length + 1]);
    CHECK_EQUAL(0, count[11 * total_length + 2]);
    CHECK_EQUAL(0, count[11 * total_length + 3]);
    // TA
    CHECK_EQUAL(0, count[12 * total_length + 0]);
    CHECK_EQUAL(0, count[12 * total_length + 1]);
    CHECK_EQUAL(0, count[12 * total_length + 2]);
    CHECK_EQUAL(0, count[12 * total_length + 3]);
    // TC
    CHECK_EQUAL(0, count[13 * total_length + 0]);
    CHECK_EQUAL(0, count[13 * total_length + 1]);
    CHECK_EQUAL(0, count[13 * total_length + 2]);
    CHECK_EQUAL(0, count[13 * total_length + 3]);
    // TG
    CHECK_EQUAL(0, count[14 * total_length + 0]);
    CHECK_EQUAL(0, count[14 * total_length + 1]);
    CHECK_EQUAL(0, count[14 * total_length + 2]);
    CHECK_EQUAL(0, count[14 * total_length + 3]);
    // TT
    CHECK_EQUAL(0, count[15 * total_length + 0]);
    CHECK_EQUAL(0, count[15 * total_length + 1]);
    CHECK_EQUAL(0, count[15 * total_length + 2]);
    CHECK_EQUAL(0, count[15 * total_length + 3]);

}


int main(int ac, char** av)
{
    return CommandLineTestRunner::RunAllTests(ac, av);
}
