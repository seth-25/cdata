//
//  main.c
//  tsgen
//
//  Created by Kostas Zoumpatianos on 3/27/12.
//  Copyright 2012 University of Trento. All rights reserved.
//
//  Modified by Karima Echihabi 15/04/2017
//  To take filename as input and write to
//  it directly instead of writing to stdout.


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm>

#include "sax/include/globals.h"
#include "sax/include/sax.h"


#define PRODUCT "TSutils - Time Series Generator\n\
Copyright (C) 2012 University of Trento\n\n"
#define STD 1   // Standard deviation





void z_normalize(float *ts, int size);

void inline z_normalize(float *ts, int size) {
    int i;
    float mean = 0;//gsl_stats_mean(ts, 1, size);
    float std = 0;//gsl_stats_sd(ts, 1, size);
    for (i=0; i<size; i++)
    {
        mean += ts[i];
    }
    mean /= size;

    for (i=0; i<size; i++)
    {
        std += (ts[i] - mean) * (ts[i] - mean);
    }
    std /= size;
    std = sqrt(std);
    for (i = 0; i < size; i++)
    {
        ts[i] = (ts[i] - mean) / std;
    }
}

float * generate (float *ts, int size, gsl_rng * r, char normalize) {
    int i;
    float x = 0, dx;

    for (i = 0; i < size; i++)
    {
        dx = gsl_ran_gaussian (r, STD); // mean=0, std=STD
        x += dx;
        ts[i] = x;
    }

    if(normalize == 1)
    {
        z_normalize(ts, size);
    }
    return ts;
}

/**
    Parses the command line arguments.
**/
void parse_args (int argc, char **argv, int *length, int *number_of_timeseries,
                 float *skew_frequency, char *normalize, char ** filename) {
    while (1)
    {
        static struct option long_options[] =  {
                {"skew-frequency", required_argument, 0, 'f'},
                {"length", required_argument, 0, 'l'},
                {"size", required_argument, 0, 's'},
                {"filename", required_argument, 0, 'o'},
                {"z-normalize", no_argument, 0, 'z'},
                {"help", no_argument, 0, 'h'},
                {NULL, 0, NULL, 0}
        };

        /* getopt_long stores the option index here. */
        int option_index = 0;

        int c = getopt_long (argc, argv, "",
                             long_options, &option_index);
        if (c == -1)
            break;
        switch (c)
        {
            case 'f':
                *skew_frequency = atof(optarg);
                break;
            case 's':
                *number_of_timeseries = atoi(optarg);
                break;
            case 'l':
                *length = atoi(optarg);
                break;
            case 'z':
                *normalize = 1;
                break;
            case 'o':
                *filename = optarg;
                break;
            case 'h':
                printf(PRODUCT);
                printf("Usage:\n\
                       \t--size XX \t\tThe number of time series to generate\n\
                       \t--length XX \t\tThe length of each time series\n\
                       \t--skew-frequency XX \tThe skewness frequency\n\
                       \t--z-normalize \t\tUse to enable z-normalization\n\
                       \t--help\n\n");
                exit(-1);
                break;
            default:
                exit(-1);
                break;
        }
    }
}

/**
    Generates a set of random time series.
**/
void generate_random_timeseries(int length, int number_of_timeseries,
                                char normalize, char * filename, bool has_saxt, char * filename1,
                                bool is_sort, bool has_ts, bool has_timestamp, int time_stamp_repeat, bool has_query_time,
                                int window_min_size, int window_max_size, int seed, char* modes) {
    // Initialize random number generation
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
//    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    gsl_rng_default_seed = (seed);
    r = gsl_rng_alloc (T);

    FILE * data_file;
    if (has_ts)
        data_file = fopen (filename, modes);


    FILE * saxt_file;
    saxt_only* saxts;
    size_t *p;
    if (has_saxt) {
        saxt_file = fopen (filename1, modes);
        saxts = (saxt_only*)malloc(sizeof(saxt_only)*number_of_timeseries);
        memset(saxts, 0, sizeof(saxt_only)*number_of_timeseries);
        p = (size_t*)malloc(sizeof(size_t)*number_of_timeseries);
        memset(p, 0, sizeof(size_t)*number_of_timeseries);
    }

    float *ts = (float*)malloc(sizeof(float) * length);

    int i;
    for (i = 1; i <= number_of_timeseries; i ++ )
    {
        generate(ts, length, r, normalize);

        if (has_ts) {
            fwrite(ts, sizeof(float), length,data_file);
            if (has_timestamp) {
                int cnt = i / time_stamp_repeat;
                if (!has_query_time) {
                    long timestamp = 1000000000 + 100 * cnt;
//                    std:: cout << "time = " << timestamp << std::endl;
                    fwrite(&timestamp, sizeof(long), 1, data_file);
                }
                else {
                    srand(seed + i);
                    int window = window_min_size + (rand() % (window_max_size - window_min_size + 1));
                    int left_cnt = cnt;
                    int right_cnt = cnt + window;
//                    std::cout << window_size << " " << left_cnt << " " << right_cnt << std::endl;
                    long l_timestamp = 1000000000 + 100 * left_cnt;
                    long r_timestamp = 1000000000 + 100 * right_cnt;
                    fwrite(&l_timestamp, sizeof(long), 1, data_file);
                    fwrite(&r_timestamp, sizeof(long), 1, data_file);
//                    std:: cout << "time = " << l_timestamp << " " << r_timestamp << std::endl;
                }
            }
        }
        if (has_saxt) {
            saxt_from_ts(ts, saxts[i].asaxt);
            p[i] = i;

        }

        if(i % 1000 == 0) {
            fprintf(stderr,"\r\x1b[m>> Generating: \x1b[36m%2.2lf%%\x1b[0m",(float) ((float)i/(float)number_of_timeseries) * 100);
        }
    }
    fprintf(stderr, "\n");

    if (has_saxt) {
        if (is_sort) {
            std::sort(saxts, saxts + number_of_timeseries + 1);
        }
        for (i = 1; i <= number_of_timeseries; i ++ ) {
            if (i % 1000 == 0)
                saxt_print(saxts[i]);
            fwrite(p+i, sizeof(size_t), 1 , saxt_file);
            fwrite(saxts[i].asaxt, sizeof(saxt_type), Bit_cardinality , saxt_file);
        }
    }


    // Finalize random number generator
    if (has_ts) fclose (data_file);
    if (has_saxt) fclose(saxt_file);
    gsl_rng_free (r);

}


//追加写入

int main(int argc, char **argv) {
//    // Parse command line arguments
//    parse_args(argc, argv, &length, &number_of_timeseries, &skew_frequency, &normalize,&filename);
//


// Initialize variables

    int length = Ts_length;                 // Length of a single time series
    int number_of_timeseries =  70000;   // Number of time series to generate    1e6=1GB
    int number_of_timeseries1 =  30000;   // Number of time series to generate    1e6=1GB

// 倾斜度
    char normalize = 1;             // Normalize or not.
    char * filename = "./query.bin";   //  File name of TS
    char * filename_saxt = "./saxt.bin";    //  File name of SAXT


    bool has_ts = 1;    //  Generate TS or not
    bool has_saxt = 0;  //  Generate SAXT or not
    bool is_sort = 0;   //  The generated SAXTs are sorted or not
    bool has_timestamp = 0; //  The generated data contains a timestamp or not
    int time_stamp_repeat = 100000; //  How many TS have the same timestamp
    bool has_query_time = 0;    // Is query or not (with two timestamps)
    int windows_max_size = 1000;    //  Randomly generated timestamp window size
    int windows_min_size = 100; //  Randomly generated timestamp window size

    int seed = 123; // The gsl seed of the previous part of the data

    int seed2 = 100;    //  The gsl seed of the latter part of the data



    fprintf(stderr,PRODUCT);
    fprintf(stderr, ">> Generating random time series...\n");
    fprintf(stderr, ">> Data Filename: %s\n", filename);
    generate_random_timeseries(length, number_of_timeseries, normalize, filename, has_saxt,
                               filename_saxt, is_sort, has_ts, has_timestamp, time_stamp_repeat,
                               has_query_time, windows_min_size, windows_max_size, seed, "w");
    generate_random_timeseries(length, number_of_timeseries1, normalize, filename, has_saxt,
                               filename_saxt, is_sort, has_ts, has_timestamp, time_stamp_repeat,
                               has_query_time, windows_min_size, windows_max_size, seed2, "a");

    fprintf(stderr, ">> Done.\n");
    return 0;
}