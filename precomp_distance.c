#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
//#include <errno.h>
#include <argp.h>

#include <hdf5_hl.h>
#include "precomp_distance_module.h"

// Prepare for argp_parse
char const *argp_program_version = "precomp_distance 0.1";
char const *argp_program_bug_address = "<example@u-tokyo.ac.jp>";
static char doc[] = "precomp_distance -- a program to compute a distance matrix among bases for SVM kernels. The distance is based on a HDF5 file containing kinetics information specified as FILE(s).";
static char args_doc[] = "FILE1 [FILE2...]";
// Keys for options without short-options
#define OPT_DATA_CONVERSION_ONLY 1
static struct argp_option options[] = {
    {0, 'k', "LENGTH", 0, "Set the length of substring to LENGTH. Default: xxx."},
    {0, 'p', "NUM", 0, "Set the position of kinetics observation from the 5' end of k-mer to NUM (1-based). Default: xxx."},
    {"chars", 'c', "STRING", 0, "Set the character set of the bases in the input kinetics file to STRING. Do not include delimiters. Default: ACGT"},
    {"threshold", 't', "VAL", 0, "Set the threshold of the rate of observed k-mers to all the possible k-mers for a distance between bases to VAL [0,1]. Default: xxx."},
    {"weight", 'w', "VAL", 0, "Set the relative weight for position-fixed complements compared with base-fixed complements to VAL [0,1]. Default: xxx."},
    {"dataset-output", 'd', "FILE", 0, "Write an IPD data set for libsvm to FILE from input data. No output by default."},
    {"output", 'o', "FILE", 0, "Write a base distance matrix to FILE. Default: standard output"},
    {"pacbio-output", 'P', "FILE", 0, "Write IPDs estimated with a PacBio model in the input data (modelPrediction) to FILE. No output by default."},
    {"data-conversion-only", OPT_DATA_CONVERSION_ONLY, 0, 0, "Only convert an input data into a libsvm data. Do not compute a k-mer distance."},
    {0}
};
struct arguments {
    char **file_paths;
    size_t file_num;
    size_t k;
    size_t p;
    char *chars;
    double count_threshold_rate;
    double position_relative_weight;
    char *dataset_output_path;
    char *distance_output_path;
    char *pacbio_model_output_path;
    int data_conversion_only;
};
// According to the manual of argp, the return type should be errno_t,
// but I couldn't use it in my environment.
static int parse_opt(int key, char *arg, struct argp_state *state){
    struct arguments *arguments = state->input;
    char *remain;
    long lparsed;
    double dparsed;
    switch(key){
        case 'k':
            lparsed = strtol(arg, &remain, 10);
            if(arg[0] == '\0' || remain[0] != '\0' || lparsed <= 0){
                fprintf(stderr, "ERROR: Invalid argument for k\n"); argp_usage(state);
            }
            arguments->k = lparsed;
            break;
        case 'p':
            lparsed = strtol(arg, &remain, 10);
            if(arg[0] == '\0' || remain[0] != '\0' || lparsed <= 0){
                fprintf(stderr, "ERROR: Invalid argument for p\n"); argp_usage(state);
            }
            arguments->p = lparsed;
            break;
        case 'c':
            arguments->chars = arg;
            break;
        case 't':
            dparsed = strtod(arg, &remain);
            if(arg[0] == '\0' || remain[0] != '\0' || dparsed < 0.0 || dparsed > 1.0){
                fprintf(stderr, "ERROR: Invalid argument for t\n"); argp_usage(state);
            }
            arguments->count_threshold_rate = dparsed;
            break;
        case 'w':
            dparsed = strtod(arg, &remain);
            if(arg[0] == '\0' || remain[0] != '\0' || dparsed < 0.0 || dparsed > 1.0){
                fprintf(stderr, "ERROR: Invalid argument for w\n"); argp_usage(state);
            }
            arguments->position_relative_weight = dparsed;
            break;
        case 'd':
            arguments->dataset_output_path = arg;
            break;
        case 'o':
            arguments->distance_output_path = arg;
            break;
        case 'P':
            arguments->pacbio_model_output_path = arg;
            break;
        case OPT_DATA_CONVERSION_ONLY:
            arguments->data_conversion_only = 1;
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

void collect_ipd_by_kmer_from_hdf5(struct arguments const *arguments, char const *file_path, size_t chars_size, size_t kmers_size, double *tMean_sum, size_t *tMean_count){
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
        if(arguments->pacbio_model_output_path != NULL){
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
        }

        // Summarize IPD
        // TODO: It may be effective to parallelize this call
        // For example, using pthread at the expence of memory usage
        collect_ipd_by_kmer(arguments->k, arguments->p, arguments->chars, tMean_buf, base_buf, (size_t)tMean_dim, tMean_sum, tMean_count, arguments->dataset_output_path, modelPrediction_buf, arguments->pacbio_model_output_path);

        // Ending process
        free(tMean_buf);
        free(base_buf[0]);
        free(base_buf);
        free(modelPrediction_buf);
        free(name);
        free(tMean_name);
        free(base_name);
    }
    H5Fclose(file_id);
    return;
}

void collect_ipd_by_kmer_from_libsvm_data(char const *file_path, double *sum, size_t *count){
    FILE *fp = fopen(file_path, "r");
    if(fp == NULL){ fprintf(stderr, "ERROR: Cannot open file: %s\n", file_path); exit(EXIT_FAILURE); }
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    char *remain;
    while((read = getline(&line, &len, fp)) != -1){
        double tMean = strtod(line, &remain);
        if(strncmp(remain, " 1:", 3) != 0){
            fprintf(stderr, "ERROR: Unexpected format in the first term: %s", line);
            exit(EXIT_FAILURE);
        }
        size_t kmer_idx = strtoull(remain + 3, &remain, 0);
        if(remain[0] != '\n'){ fprintf(stderr, "ERROR: Unexpected format in the second term: %s", line); exit(EXIT_FAILURE); }
        sum[kmer_idx] += tMean;
        ++count[kmer_idx];
    }
    free(line);
    fclose(fp);
    return;
}

int main(int argc, char **argv){
    // Default parameters
    struct arguments arguments = {
        .file_paths = NULL,
        .file_num = 0,
        .k = 12,
        .p = 9,
        .chars = "ACGT",
        .count_threshold_rate = 0.1,
        .position_relative_weight = 0.2,
        .dataset_output_path = NULL,
        .distance_output_path = NULL,
        .pacbio_model_output_path = NULL,
        .data_conversion_only = 0
    };
    // Change default parameters
    // arguments.k = 10;
    // arguments.p = 7;
    argp_parse(&argp, argc, argv, 0, 0, &arguments);
    if(arguments.dataset_output_path != NULL){
        FILE *dataset_out = fopen(arguments.dataset_output_path, "w");
        if(dataset_out == NULL){
            fprintf(stderr, "ERROR: Cannot create/truncate file: %s\n", arguments.dataset_output_path); exit(EXIT_FAILURE);
        }
        fclose(dataset_out);
    }
    fprintf(stderr, "INFO: k = %zu, p = %zu, chars = %s, threshold = %.17g, pweight = %.17g, dataset_out = %s\n", arguments.k, arguments.p, arguments.chars, arguments.count_threshold_rate, arguments.position_relative_weight, (arguments.dataset_output_path!=NULL) ? arguments.dataset_output_path : "(NONE)");
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
    size_t chars_size = strlen(arguments.chars);
    size_t kmers_size = (size_t)(pow(chars_size, arguments.k) + 0.5);
    double *tMean_sum = (double *)malloc(kmers_size * sizeof(double));
    if(tMean_sum == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for tMean_sum\n"); exit(EXIT_FAILURE); }
    memset(tMean_sum, 0, kmers_size * sizeof(double));
    size_t *tMean_count = (size_t *)malloc(kmers_size * sizeof(size_t));
    if(tMean_count == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for tMean_count\n"); exit(EXIT_FAILURE); }
    memset(tMean_count, 0, kmers_size * sizeof(size_t));

    for(size_t i = 0; i < arguments.file_num; ++i){
        size_t file_path_len = strlen(arguments.file_paths[i]);
        if(strcmp(arguments.file_paths[i] + file_path_len - 3, ".h5") == 0){
            collect_ipd_by_kmer_from_hdf5(&arguments, arguments.file_paths[i], chars_size, kmers_size, tMean_sum, tMean_count);
        }else{
            collect_ipd_by_kmer_from_libsvm_data(arguments.file_paths[i], tMean_sum, tMean_count);
        }
    }
    if(arguments.data_conversion_only == 1){ exit(EXIT_SUCCESS); }

    // Construct a distance matrix
    double *distance_matrix = (double *)malloc(arguments.k * chars_size * chars_size * sizeof(double));
    if(distance_matrix == NULL) { fprintf(stderr, "ERROR: Cannot allocate memory for distance_matrix\n"); exit(EXIT_FAILURE); }
    fill_base_distance_matrix(chars_size, arguments.k, tMean_sum, tMean_count, arguments.count_threshold_rate, arguments.position_relative_weight, distance_matrix);
    free(tMean_sum);
    free(tMean_count);

    // Write a result
    FILE *distance_out;
    if(arguments.distance_output_path == NULL){
        distance_out = stdout;
    } else {
        distance_out = fopen(arguments.distance_output_path, "w");
        if(distance_out == NULL){
            fprintf(stderr, "WARNING: Cannot create file: %s\n", arguments.distance_output_path);
            distance_out = stdout;
        }
    }
    fprintf(distance_out, "# k = %zu, p = %zu, chars = %s, threshold = %.17g, pweight = %.17g, dataset_out = %s\n", arguments.k, arguments.p, arguments.chars, arguments.count_threshold_rate, arguments.position_relative_weight, (arguments.dataset_output_path!=NULL) ? arguments.dataset_output_path : "(NONE)");
    for(size_t i = 0; i < arguments.file_num; ++i){
        fprintf(distance_out, "# file[%zu] = %s\n", i, arguments.file_paths[i]);
    }
    fprintf(distance_out, "# k, chars_size, matrix\n");
    fprintf(distance_out, "%zu\n", arguments.k);
    fprintf(distance_out, "%zu\n", chars_size);
    for(size_t position = 0; position < arguments.k; position++){
        fprintf(distance_out, "# 0-based position: %zu\n", position);
        for(size_t base1 = 0; base1 < chars_size; base1++){
            for(size_t base2 = 0; base2 < chars_size; base2++){
                size_t idx = position * chars_size * chars_size + base1 * chars_size + base2;
                fprintf(distance_out, "%.17g", distance_matrix[idx]);
                if(base2 < chars_size - 1){
                    fprintf(distance_out, "\t");
                }
            }
            fprintf(distance_out, "\n");
        }
    }
    free(distance_matrix);
    if(distance_out != stdout){ fclose(distance_out); }
    free(arguments.file_paths);
    return 0;
}
