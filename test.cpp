#include <stdint.h>
#include <cmath>
#include <CppUTest/CommandLineTestRunner.h>
#include "precomp_distance_module.h"

TEST_GROUP(kmer_ipd)
{
    // Tests for
    // int collect_ipd_by_kmer(int k, int p, char *chars, float *tMeans, char **bases, hsize_t dim, double *sum, size_t *count);
    char *chars = (char *)"ACGT";
    size_t dim = 10;
    // positive: AACGC, negative: GC0TT
    const char *bases[10]      = {"A", "T", "A", "T", "C", "\0", "G", "C", "C", "G"};
    float tMeans[10]           = {1.5, 2.2, 3.1, 4.9, 5.5, 6.3, 7.2, 8.2, 9.1, 10.3};
    float modelPredictions[10] = {1.4, 2.1, 3.0, 4.8, 5.4, 6.2, 7.1, 8.1, 9.0, 10.2};
};

TEST(kmer_ipd, k1)
{
    int k = 1;
    int p = 1;
    size_t kmers_size = (size_t)(std::pow(4, k) + 0.5);
    CHECK_EQUAL(kmers_size, 4);
    double sum[kmers_size] = {0.0};
    size_t count[kmers_size] = {0};
    collect_ipd_by_kmer(k, p, chars, tMeans, (char **)bases, dim, sum, count, NULL, NULL, NULL);
    CHECK_EQUAL((float)sum[0], 1.5f + 3.1f);
    CHECK_EQUAL((float)sum[1], 5.5f + 8.2f + 9.1f);
    CHECK_EQUAL((float)sum[2], 7.2f + 10.3f);
    CHECK_EQUAL((float)sum[3], 2.2f + 4.9f);
    CHECK_EQUAL(count[0], 2);
    CHECK_EQUAL(count[1], 3);
    CHECK_EQUAL(count[2], 2);
    CHECK_EQUAL(count[3], 2);
}

TEST(kmer_ipd, k2p1)
{
    int k = 2;
    int p = 1;
    size_t kmers_size = (size_t)(std::pow(4, k) + 0.5);
    CHECK_EQUAL(kmers_size, 16);
    double sum[kmers_size] = {0.0};
    size_t count[kmers_size] = {0};
    collect_ipd_by_kmer(k, p, chars, tMeans, (char **)bases, dim, sum, count, NULL, NULL, NULL);
    CHECK_EQUAL((float)sum[0], 1.5f); // AA
    CHECK_EQUAL(count[0], 1);
    CHECK_EQUAL((float)sum[1], 3.1f); // AC
    CHECK_EQUAL(count[1], 1);
    CHECK_EQUAL((float)sum[2], 0.0f); // AG
    CHECK_EQUAL(count[2], 0);
    CHECK_EQUAL((float)sum[3], 0.0f); // AT
    CHECK_EQUAL(count[3], 0);
    CHECK_EQUAL((float)sum[6], 5.5f); // CG
    CHECK_EQUAL(count[6], 1);
    CHECK_EQUAL((float)sum[9], 7.2f + 10.3f); // GC
    CHECK_EQUAL(count[9], 2);
    CHECK_EQUAL((float)sum[15], 4.9f); // TT
    CHECK_EQUAL(count[15], 1);
}

TEST(kmer_ipd, k2p2)
{
    int k = 2;
    int p = 2;
    size_t kmers_size = (size_t)(std::pow(4, k) + 0.5);
    CHECK_EQUAL(kmers_size, 16);
    double sum[kmers_size] = {0.0};
    size_t count[kmers_size] = {0};
    collect_ipd_by_kmer(k, p, chars, tMeans, (char **)bases, dim, sum, count, "test.tmp.data_for_libsvm", modelPredictions, "test.tmp.pacbio_data");
    CHECK_EQUAL((float)sum[0], 3.1f); // AA
    CHECK_EQUAL(count[0], 1);
    CHECK_EQUAL((float)sum[1], 5.5f); // AC
    CHECK_EQUAL(count[1], 1);
    CHECK_EQUAL((float)sum[2], 0.0f); // AG
    CHECK_EQUAL(count[2], 0);
    CHECK_EQUAL((float)sum[3], 0.0f); // AT
    CHECK_EQUAL(count[3], 0);
    CHECK_EQUAL((float)sum[6], 7.2f); // CG
    CHECK_EQUAL(count[6], 1);
    CHECK_EQUAL((float)sum[9], 9.1f + 8.2f); // GC
    CHECK_EQUAL(count[9], 2);
    CHECK_EQUAL((float)sum[15], 2.2f); // TT
    CHECK_EQUAL(count[15], 1);
}

TEST(kmer_ipd, k4p3)
{
    int k = 4;
    int p = 3;
    size_t kmers_size = (size_t)(std::pow(4, k) + 0.5);
    CHECK_EQUAL(kmers_size, 256);
    double sum[kmers_size] = {0.0};
    size_t count[kmers_size] = {0};
    collect_ipd_by_kmer(k, p, chars, tMeans, (char **)bases, dim, sum, count, NULL, NULL, NULL);
    CHECK_EQUAL((float)sum[0], 0.0f); // AAAA
    CHECK_EQUAL(count[0], 0);
    CHECK_EQUAL((float)sum[5], 0.0f); // AACC
    CHECK_EQUAL(count[5], 0);
    CHECK_EQUAL((float)sum[6], 5.5f); // AACG
    CHECK_EQUAL(count[6], 1);
    CHECK_EQUAL((float)sum[7], 0.0f); // AACT
    CHECK_EQUAL(count[7], 0);
    CHECK_EQUAL((float)sum[24], 0.0f); // ACGA
    CHECK_EQUAL(count[24], 0);
    CHECK_EQUAL((float)sum[25], 7.2f); // ACGC
    CHECK_EQUAL(count[25], 1);
    CHECK_EQUAL((float)sum[26], 0.0f); // ACGG
    CHECK_EQUAL(count[26], 0);
    CHECK_EQUAL((float)sum[241], 0.0f); // TTAC
    CHECK_EQUAL(count[241], 0);
    CHECK_EQUAL((float)sum[255], 0.0f); // TTTT
    CHECK_EQUAL(count[255], 0);
}

TEST(kmer_ipd, k4p3update)
{
    int k = 4;
    int p = 3;
    size_t kmers_size = (size_t)(std::pow(4, k) + 0.5);
    CHECK_EQUAL(kmers_size, 256);
    double sum[kmers_size] = {0.0};
    sum[24] = 240.8;
    sum[25] = 250.9;
    size_t count[kmers_size] = {0};
    count[24] = 4200;
    count[25] = 25005;
    collect_ipd_by_kmer(k, p, chars, tMeans, (char **)bases, dim, sum, count, NULL, NULL, NULL);
    CHECK_EQUAL((float)sum[0], 0.0f); // AAAA
    CHECK_EQUAL(count[0], 0);
    CHECK_EQUAL((float)sum[5], 0.0f); // AACC
    CHECK_EQUAL(count[5], 0);
    CHECK_EQUAL((float)sum[6], 5.5f); // AACG
    CHECK_EQUAL(count[6], 1);
    CHECK_EQUAL((float)sum[7], 0.0f); // AACT
    CHECK_EQUAL(count[7], 0);
    CHECK_EQUAL(sum[24], 240.8); // ACGA
    CHECK_EQUAL(count[24], 4200);
    CHECK_EQUAL(sum[25], 250.9 + 7.2f); // ACGC
    CHECK_EQUAL(count[25], 25006);
    CHECK_EQUAL((float)sum[26], 0.0f); // ACGG
    CHECK_EQUAL(count[26], 0);
    CHECK_EQUAL((float)sum[241], 0.0f); // TTAC
    CHECK_EQUAL(count[241], 0);
    CHECK_EQUAL((float)sum[255], 0.0f); // TTTT
    CHECK_EQUAL(count[255], 0);
}

TEST_GROUP(insert_base){
    // size_t insert_base_into_context(size_t context, size_t chars_size, size_t position, size_t base);
};

TEST(insert_base, test1){
    // Assuming chars == "ACGT"
    size_t context = 241; // TTAC
    size_t chars_size = 4;
    int k = 5;
    size_t position = 1; // 0-based index
    size_t base = 2; // G
    size_t ret = insert_base_into_context(context, chars_size, k, position, base);
    CHECK_EQUAL(ret, 945); // TGTAC
}

TEST(insert_base, test2){
    // Assuming chars == "ACGT"
    size_t context = 241; // TTAC
    size_t chars_size = 4;
    int k = 5;
    size_t position = 4; // 0-based index
    size_t base = 2; // G
    size_t ret = insert_base_into_context(context, chars_size, k, position, base);
    CHECK_EQUAL(ret, 966); // TTACG
}

TEST(insert_base, test3){
    // Assuming chars == "ACGT"
    size_t context = 241; // TTAC
    size_t chars_size = 4;
    int k = 5;
    size_t position = 0; // 0-based index
    size_t base = 2; // G
    size_t ret = insert_base_into_context(context, chars_size, k, position, base);
    CHECK_EQUAL(ret, 753); // GTTAC
}

TEST_GROUP(compute_distance){
    // Tests for
    // double compute_base_distance(size_t position, size_t base1, size_t base2, size_t chars_size, int k, double *ipd_sum, size_t *count, double count_threshold_rate);
    size_t chars_size = 4;
    int k = 4;
    double ipd_sum[256] = {0.0};
    size_t count[256] = {0};
    double p1answer;
    double p3answer;
    size_t min_valid_context_count;

    void setup(){
        min_valid_context_count = SIZE_MAX;
        ipd_sum[68] = 17.17; // CACA
        count[68] = 6;
        ipd_sum[69] = 20.2; // CACC
        count[69] = 7;
        ipd_sum[84] = 14.11; // CCCA
        count[84] = 12;
        ipd_sum[85] = 10.3; // CCCC
        count[85] = 11;
        ipd_sum[89] = 8.2; // CCGC
        count[89] = 10;
        p1answer = std::sqrt((((ipd_sum[68]/count[68] - ipd_sum[84]/count[84]) * (ipd_sum[68]/count[68] - ipd_sum[84]/count[84])) + ((ipd_sum[69]/count[69] - ipd_sum[85]/count[85]) * (ipd_sum[69]/count[69] - ipd_sum[85]/count[85]))) / 2);
        p3answer = std::sqrt((((ipd_sum[68]/count[68] - ipd_sum[69]/count[69]) * (ipd_sum[68]/count[68] - ipd_sum[69]/count[69])) + ((ipd_sum[84]/count[84] - ipd_sum[85]/count[85]) * (ipd_sum[84]/count[84] - ipd_sum[85]/count[85]))) / 2);
    }
};

TEST(compute_distance, p1pass){
    size_t position = 1;
    size_t base1 = 0; // A
    size_t base2 = 1; // C
    double count_threshold_rate = 0.031; // Actual rate should be 0.0312
    double ret = compute_base_distance(position, base1, base2, chars_size, k, ipd_sum, count, count_threshold_rate, &min_valid_context_count);
    CHECK_EQUAL(p1answer, ret);
    CHECK_EQUAL(2, min_valid_context_count);
}

TEST(compute_distance, p1fail){
    size_t position = 1;
    size_t base1 = 0; // A
    size_t base2 = 1; // C
    double count_threshold_rate = 0.032; // Actual rate should be 0.0312
    double ret = compute_base_distance(position, base1, base2, chars_size, k, ipd_sum, count, count_threshold_rate, &min_valid_context_count);
    CHECK_EQUAL(-1.0, ret);
    CHECK_EQUAL(2, min_valid_context_count);
}

TEST(compute_distance, p3pass){
    size_t position = 3;
    size_t base1 = 0; // A
    size_t base2 = 1; // C
    double count_threshold_rate = 0.031; // Actual rate should be 0.0312
    double ret = compute_base_distance(position, base1, base2, chars_size, k, ipd_sum, count, count_threshold_rate, &min_valid_context_count);
    CHECK_EQUAL(p3answer, ret);
    CHECK_EQUAL(2, min_valid_context_count);
}

TEST(compute_distance, fill_distance_matrix){
    // Test for
    // int fill_base_distance_matrix(size_t chars_size, int k, double *ipd_sum, size_t *count, double count_threshold_rate, double position_relative_weight, double *matrix);
    double matrix[256] = {0.0};
    double count_threshold_rate = 0.031; // A value between 1/64 and 2/64
    double position_relative_weight = 0.2;
    fill_base_distance_matrix(chars_size, k, ipd_sum, count, count_threshold_rate, position_relative_weight, matrix);
    CHECK_EQUAL(matrix[3*chars_size*chars_size + 2*chars_size + 3], -HUGE_VAL); // pos:3, base1:G, base2:T
    CHECK_EQUAL(matrix[3*chars_size*chars_size + 3*chars_size + 3], 0.0); // pos:3, base1:T, base2:T
    CHECK_EQUAL(matrix[1*chars_size*chars_size + 0*chars_size + 1], p1answer); // pos:1, base1:A, base2:C
    CHECK_EQUAL(matrix[1*chars_size*chars_size + 1*chars_size + 0], p1answer); // pos:1, base1:C, base2:A
    CHECK_EQUAL(matrix[2*chars_size*chars_size + 0*chars_size + 1], -HUGE_VAL); // pos:2, base1:A, base2:C
    CHECK_EQUAL(matrix[2*chars_size*chars_size + 1*chars_size + 2], -HUGE_VAL); // pos:2, base1:C, base2:G
    // Check using another threshold
    double matrix2[256] = {0.0};
    count_threshold_rate = 0.01; // less than 1/64
    fill_base_distance_matrix(chars_size, k, ipd_sum, count, count_threshold_rate, position_relative_weight, matrix2);
    CHECK(matrix2[2*chars_size*chars_size + 1*chars_size + 2] > 0.0); // pos:2, base1:C, base2:G
    double p2ac_complement = position_relative_weight * matrix2[2*chars_size*chars_size + 1*chars_size + 2] + (1-position_relative_weight) * (matrix2[1*chars_size*chars_size + 1] + matrix2[3*chars_size*chars_size + 1*chars_size]) / 2;
    CHECK_EQUAL(matrix2[2*chars_size*chars_size + 0*chars_size + 1], p2ac_complement); // pos:2, base1:A, base2:C
}

int main(int ac, char** av)
{
    return CommandLineTestRunner::RunAllTests(ac, av);
}

