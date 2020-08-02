#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
//#include <errno.h>
#include <argp.h>

#include <hdf5_hl.h>
#include "collect_ipd_module.h"

// Prepare for argp_parse
char const *argp_program_version = "collect_ipd 1.0";
char const *argp_program_bug_address = "<example@u-tokyo.ac.jp>";
static char doc[] = "collect_ipd -- a program to collect IPDs and model prediction in HDF5 files specified as FILE(s).";
static char args_doc[] = "FILE1 [FILE2...]";
// Keys for options without short-options
#define OPT_DATA_CONVERSION_ONLY 1
static struct argp_option options[] = {
    {0, 'k', "LENGTH", 0, "Set the length of substring (k-mer) to LENGTH. Default: 2."},
    {0, 'l', "LENGTH", 0, "Set the outside length of k-mers. Default: 20."},
    {"chars", 'c', "STRING", 0, "Set the character set of the bases in the input kinetics file to STRING. Do not include delimiters. Default: ACGT"},
    {"threshold", 't', "INTEGER", 0, "Set the threshold of coverage of observed k-mers. Default: 25."},
    {"output", 'o', "FILE", 0, "Write IPD sum per k-mer to FILE. Default: standard output"},
    {0}
};
struct arguments {
    char **file_paths;
    size_t file_num;
    size_t k;
    size_t outside_length;
    char *chars;
    size_t coverage_threshold;
    char *output_path;
};
// According to the manual of argp, the return type should be errno_t,
// but I couldn't use it in my environment.
static int parse_opt(int key, char *arg, struct argp_state *state){
    struct arguments *arguments = state->input;
    char *remain;
    long lparsed;
    switch(key){
        case 'k':
            lparsed = strtol(arg, &remain, 10);
            if(arg[0] == '\0' || remain[0] != '\0' || lparsed <= 0){
                fprintf(stderr, "ERROR: Invalid argument for k\n"); argp_usage(state);
            }
            arguments->k = lparsed;
            break;
        case 'l':
            lparsed = strtol(arg, &remain, 10);
            if(arg[0] == '\0' || remain[0] != '\0' || lparsed < 0){
                fprintf(stderr, "ERROR: Invalid argument for l\n"); argp_usage(state);
            }
            arguments->outside_length = lparsed;
            break;
        case 'c':
            arguments->chars = arg;
            break;
        case 't':
            lparsed = strtol(arg, &remain, 10);
            if(arg[0] == '\0' || remain[0] != '\0' || lparsed < 0){
                fprintf(stderr, "ERROR: Invalid argument for t\n"); argp_usage(state);
            }
            arguments->coverage_threshold = lparsed;
            break;
        case 'o':
            arguments->output_path = arg;
            break;
        case ARGP_KEY_ARG:
            arguments->file_num++;
            arguments->file_paths = (char **)realloc(arguments->file_paths, sizeof(char *) * arguments->file_num);
            if(arguments->file_paths == NULL){ fprintf(stderr, "ERROR: Cannot realloc file_paths\n"); }
            arguments->file_paths[arguments->file_num - 1] = arg;
            break;
        case ARGP_KEY_END:
            if(state->arg_num < 1){
                fprintf(stderr, "ERROR: Too few arguments\n"); argp_usage(state);
            }
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}
static struct argp argp = {options, parse_opt, args_doc, doc};

// Write IPD data per k-mer
// Column: k-mer index, k-mer string, position (1 == start of k-mer), chromosome name, IPD sum, squared IPD sum, model prediction sum, squared model prediction sum, count
void write_ipd_by_kmer(size_t const k, size_t const outside_length, size_t const chars_size, char const *chars, char const *chromosome_name, size_t const file_idx,
        double const *tMean_sum, double const *tMean_sq_sum, double const *tMean_log2_sum, double const *tMean_log2_sq_sum,
        double const *prediction_sum, double const *prediction_sq_sum, double const *prediction_log2_sum, double const *prediction_log2_sq_sum, size_t const *count,
        int const print_header, FILE *output) {
    size_t kmers_size = (size_t)(pow(chars_size, k) + 0.5);
    size_t total_length = k + 2 * outside_length;
    char *kmer_string = (char *)malloc((k + 1) * sizeof(char));
    if(kmer_string == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for kmer_string\n"); exit(EXIT_FAILURE); }
    kmer_string[k] = '\0';
    if(print_header == 1) {
        fprintf(output, "kmer_string,kmer_number,position,chromosome,file_index,ipd_sum,ipd_sq_sum,log2_ipd_sum,log2_ipd_sq_sum,prediction_sum,prediction_sq_sum,log2_prediction_sum,log2_prediction_sq_sum,count\n");
    }
    for (size_t kmer = 0; kmer < kmers_size; ++kmer) {
        size_t kmer_tmp = kmer;
        for (int i = (int)k - 1; i >= 0; --i) {
            kmer_string[i] = chars[kmer_tmp % chars_size];
            kmer_tmp /= chars_size;
        }
        for (size_t i = 0; i < total_length; ++i) {
            size_t idx = kmer * total_length + i;
            fprintf(output, "%s,%zu,%d,%s,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%zu\n",
                    kmer_string, kmer, (int)i - (int)outside_length + 1, chromosome_name, file_idx, tMean_sum[idx], tMean_sq_sum[idx], tMean_log2_sum[idx], tMean_log2_sq_sum[idx],
                    prediction_sum[idx], prediction_sq_sum[idx], prediction_log2_sum[idx], prediction_log2_sq_sum[idx], count[idx]);
        }
    }
    free(kmer_string);
    return;
}

void collect_ipd_by_kmer_from_hdf5(char const *file_path, size_t const file_index, size_t const k, size_t const outside_length,
        size_t const chars_size, char const *chars, size_t const kmers_size, size_t const coverage_threshold, FILE *output,
        double *tMean_sum, double *tMean_sq_sum, double *tMean_log2_sum, double *tMean_log2_sq_sum,
        double *prediction_sum, double *prediction_sq_sum, double *prediction_log2_sum, double *prediction_log2_sq_sum, size_t *count){
    if(sizeof(hsize_t) < sizeof(size_t)){
        fprintf(stderr, "WARNING: sizeof(hsize_t) == %zu < sizeof(size_t) == %zu: the result may be incorrect\n", sizeof(hsize_t), sizeof(size_t));
    }
    herr_t hstatus;
    hid_t file_id = H5Fopen(file_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if(file_id < 0) { fprintf(stderr, "ERROR: Cannot open file in HDF5 format: %s\n", file_path); exit(EXIT_FAILURE); }
    H5G_info_t ginfo;
    H5Gget_info_by_name(file_id, "/", &ginfo, H5P_DEFAULT);
    fprintf(stderr, "INFO: # chromosomes: %d\n", (int)ginfo.nlinks);
    // Scan data sets for each chromosome
    for (size_t i = 0; i < ginfo.nlinks; i++){
        ssize_t name_size = H5Lget_name_by_idx(file_id, "/", H5_INDEX_NAME, H5_ITER_NATIVE, i, NULL, 0, H5P_DEFAULT);
        //printf("%zd\n", name_size);
        if(name_size < 0) { fprintf(stderr, "ERROR: Cannot get name by idx: %zu\n", i); exit(EXIT_FAILURE); }
        if(name_size >= SIZE_MAX / 2) { fprintf(stderr, "Name is too long: %zd\n", name_size); exit(EXIT_FAILURE); }
        // Note: sizeof(char) == 1
        char *name = (char *)malloc(name_size + 1);
        if(name == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for name\n"); exit(EXIT_FAILURE); }
        H5Lget_name_by_idx(file_id, "/", H5_INDEX_NAME, H5_ITER_NATIVE, i, name, name_size + 1, H5P_DEFAULT);
        //printf("%s\n", name);
        // Make a dataset name to access the dataset
        // dataset: tMean
        char const *tMean = "tMean";
        char *tMean_name = (char *)malloc(name_size + strlen(tMean) + 3);
        if(tMean_name == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for tMean_name\n"); exit(EXIT_FAILURE); }
        sprintf(tMean_name, "/%s/%s", name, tMean);
        int rank = 0;
        H5LTget_dataset_ndims(file_id,tMean_name,&rank);
        if(rank != 1) { fprintf(stderr, "ERROR: Rank is not 1; observed: %d, path: %s\n", rank, tMean_name); exit(EXIT_FAILURE); }
        hid_t tMean_dset_id = H5Dopen(file_id, tMean_name, H5P_DEFAULT);
        if(tMean_dset_id < 0) { fprintf(stderr, "ERROR: Failure in opening %s\n", tMean_name); exit(EXIT_FAILURE); }
        hid_t tMean_dtype_id = H5Dget_type(tMean_dset_id);
        if(H5Tequal(tMean_dtype_id, H5T_NATIVE_FLOAT) <= 0) { fprintf(stderr, "ERROR: Dataset tMean is not H5T_NATIVE_FLOAT type. Check using h5ls.\n"); exit(EXIT_FAILURE); }
        H5Tclose(tMean_dtype_id);
        H5Dclose(tMean_dset_id);
        hsize_t tMean_dim = 0;
        H5LTget_dataset_info(file_id,tMean_name,&tMean_dim,NULL,NULL);
        fprintf(stderr, "INFO: chromosome: %s, length: %llu\n", name, tMean_dim);
        float *tMean_buf = (float *)malloc(sizeof(float) * tMean_dim);
        if(tMean_buf == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for tMean_buf\n"); exit(EXIT_FAILURE); }
        hstatus = H5LTread_dataset_float(file_id, tMean_name, tMean_buf);
        if(hstatus < 0) { fprintf(stderr, "ERROR: Failure in reading Dataset %s\n", tMean_name); exit(EXIT_FAILURE); }
        //if(strcmp(name, "chrIV") == 0) printf("ret:%d, val:%f\n", ret, buf[28173040 - 2]);
        // dataset: base
        char const *base = "base";
        char *base_name = (char *)malloc(name_size + strlen(base) + 3);
        if(base_name == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for base_name\n"); exit(EXIT_FAILURE); }
        sprintf(base_name, "/%s/%s", name, base);
        H5LTget_dataset_ndims(file_id,base_name,&rank);
        if(rank != 1) { fprintf(stderr, "ERROR: Rank is not 1; observed: %d, path: %s\n", rank, base_name); exit(EXIT_FAILURE); }
        hid_t base_dset_id = H5Dopen(file_id, base_name, H5P_DEFAULT);
        if(base_dset_id < 0) { fprintf(stderr, "ERROR: Failure in opening %s\n", base_name); exit(EXIT_FAILURE); }
        hid_t base_dtype_id = H5Dget_type(base_dset_id);
        if(H5Tget_class(base_dtype_id) != H5T_STRING) { fprintf(stderr, "ERROR: Dataset base is not H5T_STRING class. Check the input file or PacBio specification.\n"); exit(EXIT_FAILURE); }
        size_t base_len = H5Tget_size(base_dtype_id);
        H5Tclose(base_dtype_id);
        if(base_len != 1) { fprintf(stderr, "ERROR: Length of base string is not 1; observed: %zu (Dataset: %s)\n", base_len, base_name); exit(EXIT_FAILURE); }
        // Make room for null terminator. Now, base_len == 2
        base_len++;
        hsize_t base_dim = 0;
        H5LTget_dataset_info(file_id,tMean_name,&base_dim,NULL,NULL);
        if(tMean_dim != base_dim) { fprintf(stderr, "ERROR: Dataset dimension is inconsistent\n"); exit(EXIT_FAILURE); }
        char **base_buf = (char **)malloc(base_dim * sizeof(char *));
        if(base_buf == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for base_buf\n"); exit(EXIT_FAILURE); }
        base_buf[0] = (char *)malloc(base_dim * base_len);
        if(base_buf[0] == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for base_buf[0]\n"); exit(EXIT_FAILURE); }
        for(size_t j = 1; j < base_dim; j++) { base_buf[j] = base_buf[0] + j * base_len; }
        hid_t base_memtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(base_memtype, base_len);
        hstatus = H5Dread(base_dset_id, base_memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, base_buf[0]);
        //hstatus = H5LTread_dataset_string(file_id, base_name, base_buf[0]);
        if(hstatus < 0) { fprintf(stderr, "ERROR: Failure in reading Dataset %s\n", base_name); exit(EXIT_FAILURE); }
        //if(strcmp(name, "chrIV") == 0) printf("base:%.100s\n", base_buf[28173040 - 2]);
        H5Tclose(base_memtype);
        H5Dclose(base_dset_id);

        // dataset: modelPrediction
        float *modelPrediction_buf = NULL;
        char const *modelPrediction = "modelPrediction";
        char *modelPrediction_name = (char *)malloc(name_size + strlen(modelPrediction) + 3);
        if(modelPrediction_name == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for modelPrediction_name\n"); exit(EXIT_FAILURE); }
        sprintf(modelPrediction_name, "/%s/%s", name, modelPrediction);
        H5LTget_dataset_ndims(file_id,modelPrediction_name,&rank);
        if(rank != 1) { fprintf(stderr, "ERROR: Rank is not 1; observed: %d, path: %s\n", rank, modelPrediction_name); exit(EXIT_FAILURE); }
        hid_t modelPrediction_dset_id = H5Dopen(file_id, modelPrediction_name, H5P_DEFAULT);
        if(modelPrediction_dset_id < 0) { fprintf(stderr, "ERROR: Failure in opening %s\n", modelPrediction_name); exit(EXIT_FAILURE); }
        hid_t modelPrediction_dtype_id = H5Dget_type(modelPrediction_dset_id);
        if(H5Tequal(modelPrediction_dtype_id, H5T_NATIVE_FLOAT) <= 0) { fprintf(stderr, "ERROR: Dataset modelPrediction is not H5T_NATIVE_FLOAT type. Check using h5ls.\n"); exit(EXIT_FAILURE); }
        H5Tclose(modelPrediction_dtype_id);
        H5Dclose(modelPrediction_dset_id);
        hsize_t modelPrediction_dim = 0;
        H5LTget_dataset_info(file_id,modelPrediction_name,&modelPrediction_dim,NULL,NULL);
        if(tMean_dim != modelPrediction_dim) { fprintf(stderr, "ERROR: Dataset dimension is inconsistent\n"); exit(EXIT_FAILURE); }
        modelPrediction_buf = (float *)malloc(sizeof(float) * modelPrediction_dim);
        if(modelPrediction_buf == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for modelPrediction_buf\n"); exit(EXIT_FAILURE); }
        hstatus = H5LTread_dataset_float(file_id, modelPrediction_name, modelPrediction_buf);
        if(hstatus < 0) { fprintf(stderr, "ERROR: Failure in reading Dataset %s\n", modelPrediction_name); exit(EXIT_FAILURE); }
        free(modelPrediction_name);

        // dataset: coverage
        unsigned int *coverage_buf = NULL;
        char const *coverage = "coverage";
        char *coverage_name = (char *)malloc(name_size + strlen(coverage) + 3);
        if(coverage_name == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for coverage_name\n"); exit(EXIT_FAILURE); }
        sprintf(coverage_name, "/%s/%s", name, coverage);
        H5LTget_dataset_ndims(file_id,coverage_name,&rank);
        if(rank != 1) { fprintf(stderr, "ERROR: Rank is not 1; observed: %d, path: %s\n", rank, coverage_name); exit(EXIT_FAILURE); }
        hid_t coverage_dset_id = H5Dopen(file_id, coverage_name, H5P_DEFAULT);
        if(coverage_dset_id < 0) { fprintf(stderr, "ERROR: Failure in opening %s\n", coverage_name); exit(EXIT_FAILURE); }
        hid_t coverage_dtype_id = H5Dget_type(coverage_dset_id);
        if(H5Tequal(coverage_dtype_id, H5T_NATIVE_UINT) <= 0) { fprintf(stderr, "ERROR: Dataset coverage is not H5T_NATIVE_UINT type. Check using h5ls.\n"); exit(EXIT_FAILURE); }
        H5Tclose(coverage_dtype_id);
        H5Dclose(coverage_dset_id);
        hsize_t coverage_dim = 0;
        H5LTget_dataset_info(file_id,coverage_name,&coverage_dim,NULL,NULL);
        if(tMean_dim != coverage_dim) { fprintf(stderr, "ERROR: Dataset dimension is inconsistent\n"); exit(EXIT_FAILURE); }
        coverage_buf = (unsigned int *)malloc(sizeof(unsigned int) * coverage_dim);
        if(coverage_buf == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for coverage_buf\n"); exit(EXIT_FAILURE); }
        hstatus = H5LTread_dataset(file_id, coverage_name, H5T_NATIVE_UINT, coverage_buf);
        if(hstatus < 0) { fprintf(stderr, "ERROR: Failure in reading Dataset %s\n", coverage_name); exit(EXIT_FAILURE); }
        free(coverage_name);

        // Summarize IPD
        // TODO: It may be effective to parallelize this call
        // For example, using pthread at the expence of memory usage
        int check_outside_coverage = 1;
        collect_ipd_by_kmer(k, chars, tMean_buf, base_buf, (size_t)tMean_dim, tMean_sum, tMean_sq_sum, tMean_log2_sum, tMean_log2_sq_sum,
                prediction_sum, prediction_sq_sum, prediction_log2_sum, prediction_log2_sq_sum, count, modelPrediction_buf, coverage_buf, coverage_threshold, outside_length, check_outside_coverage);

        // Write data per chromosome
        int print_header = (i == 0) ? 1 : 0;
        write_ipd_by_kmer(k, outside_length, chars_size, chars, name, file_index,
                tMean_sum, tMean_sq_sum, tMean_log2_sum, tMean_log2_sq_sum, prediction_sum, prediction_sq_sum, prediction_log2_sum, prediction_log2_sq_sum, count, print_header, output);

        // Ending process
        free(tMean_buf);
        free(base_buf[0]);
        free(base_buf);
        free(modelPrediction_buf);
        free(coverage_buf);
        free(name);
        free(tMean_name);
        free(base_name);

        // Reset tMean_sum and so on
        size_t total_length = k + 2 * outside_length;
        memset(tMean_sum, 0, total_length * sizeof(double));
        memset(tMean_sq_sum, 0, total_length * sizeof(double));
        memset(tMean_log2_sum, 0, total_length * sizeof(double));
        memset(tMean_log2_sq_sum, 0, total_length * sizeof(double));
        memset(prediction_sum, 0, total_length * sizeof(double));
        memset(prediction_sq_sum, 0, total_length * sizeof(double));
        memset(prediction_log2_sum, 0, total_length * sizeof(double));
        memset(prediction_log2_sq_sum, 0, total_length * sizeof(double));
        memset(count, 0, total_length * sizeof(size_t));
    }
    H5Fclose(file_id);
    return;
}

int main(int argc, char **argv){
    // Default parameters
    struct arguments arguments = {
        .file_paths = NULL,
        .file_num = 0,
        .k = 2,
        .outside_length = 20,
        .chars = "ACGT",
        .coverage_threshold = 25,
        .output_path = NULL,
    };
    // Change default parameters
    // arguments.k = 10;
    argp_parse(&argp, argc, argv, 0, 0, &arguments);
    fprintf(stderr, "INFO: k = %zu, outside_length = %zu, chars = %s, coverage_threshold = %zu, output_path = %s\n",
            arguments.k, arguments.outside_length, arguments.chars, arguments.coverage_threshold, (arguments.output_path!=NULL) ? arguments.output_path : "(NONE)");
    for(size_t i = 0; i < arguments.file_num; ++i){
        fprintf(stderr, "INFO: file[%zu] = %s\n", i, arguments.file_paths[i]);
        FILE *tmp_fp = fopen(arguments.file_paths[i], "r");
        if(tmp_fp == NULL){
            fprintf(stderr, "ERROR: Cannot open file: %s\n", arguments.file_paths[i]);
            exit(EXIT_FAILURE);
        }else{
            fclose(tmp_fp);
        }
    }
    FILE *output;
    if(arguments.output_path == NULL){
        output = stdout;
    } else {
        output = fopen(arguments.output_path, "w");
        if(output == NULL){
            fprintf(stderr, "ERROR: Cannot create/truncate file: %s\n", arguments.output_path); exit(EXIT_FAILURE);
        }
    }
    size_t chars_size = strlen(arguments.chars);
    size_t kmers_size = (size_t)(pow(chars_size, arguments.k) + 0.5);
    size_t total_length = kmers_size * (arguments.k + 2 * arguments.outside_length);
    double *tMean_sum = (double *)malloc(total_length * sizeof(double));
    if(tMean_sum == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for tMean_sum\n"); exit(EXIT_FAILURE); }
    memset(tMean_sum, 0, total_length * sizeof(double));
    double *tMean_sq_sum = (double *)malloc(total_length * sizeof(double));
    if(tMean_sq_sum == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for tMean_sq_sum\n"); exit(EXIT_FAILURE); }
    memset(tMean_sq_sum, 0, total_length * sizeof(double));
    double *tMean_log2_sum = (double *)malloc(total_length * sizeof(double));
    if(tMean_log2_sum == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for tMean_log2_sum\n"); exit(EXIT_FAILURE); }
    memset(tMean_log2_sum, 0, total_length * sizeof(double));
    double *tMean_log2_sq_sum = (double *)malloc(total_length * sizeof(double));
    if(tMean_log2_sq_sum == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for tMean_log2_sq_sum\n"); exit(EXIT_FAILURE); }
    memset(tMean_log2_sq_sum, 0, total_length * sizeof(double));
    double *prediction_sum = (double *)malloc(total_length * sizeof(double));
    if(prediction_sum == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for prediction_sum\n"); exit(EXIT_FAILURE); }
    memset(prediction_sum, 0, total_length * sizeof(double));
    double *prediction_sq_sum = (double *)malloc(total_length * sizeof(double));
    if(prediction_sq_sum == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for prediction_sq_sum\n"); exit(EXIT_FAILURE); }
    memset(prediction_sq_sum, 0, total_length * sizeof(double));
    double *prediction_log2_sum = (double *)malloc(total_length * sizeof(double));
    if(prediction_log2_sum == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for prediction_log2_sum\n"); exit(EXIT_FAILURE); }
    memset(prediction_log2_sum, 0, total_length * sizeof(double));
    double *prediction_log2_sq_sum = (double *)malloc(total_length * sizeof(double));
    if(prediction_log2_sq_sum == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for prediction_log2_sq_sum\n"); exit(EXIT_FAILURE); }
    memset(prediction_log2_sq_sum, 0, total_length * sizeof(double));
    size_t *count = (size_t *)malloc(total_length * sizeof(size_t));
    if(count == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for count\n"); exit(EXIT_FAILURE); }
    memset(count, 0, total_length * sizeof(size_t));

    for(size_t i = 0; i < arguments.file_num; ++i){
        size_t file_path_len = strlen(arguments.file_paths[i]);
        if(strcmp(arguments.file_paths[i] + file_path_len - 3, ".h5") != 0){
            fprintf(stderr, "WARNING: %s may not be a HDF5 file. Continuing.", arguments.file_paths[i]);
        }
        collect_ipd_by_kmer_from_hdf5(arguments.file_paths[i], i, arguments.k, arguments.outside_length, chars_size, arguments.chars, kmers_size, arguments.coverage_threshold, output,
                tMean_sum, tMean_sq_sum, tMean_log2_sum, tMean_log2_sq_sum, prediction_sum, prediction_sq_sum, prediction_log2_sum, prediction_log2_sq_sum, count);
    }

    fclose(output);
    free(arguments.file_paths);
    free(tMean_sum);
    free(tMean_sq_sum);
    free(tMean_log2_sum);
    free(tMean_log2_sq_sum);
    free(prediction_sum);
    free(prediction_sq_sum);
    free(prediction_log2_sum);
    free(prediction_log2_sq_sum);
    free(count);
    return 0;
}
