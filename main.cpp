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
#include <readline/readline.h>
#include <getopt.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <zconf.h>
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
                                char normalize, int repetition, char * filename, bool issaxt,char * filename1,
                                bool issort, bool ists, bool hastimestamp) {
    // Initialize random number generation
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r = gsl_rng_alloc (T);
    FILE * data_file;
    if (ists)
    data_file = fopen (filename,"w");


    FILE * saxt_file;
    saxt_only* saxts;
    size_t *p;
    if (issaxt) {
        saxt_file = fopen (filename1,"w");
        saxts = (saxt_only*)malloc(sizeof(saxt_only)*number_of_timeseries);
        memset(saxts, 0, sizeof(saxt_only)*number_of_timeseries);
        p = (size_t*)malloc(sizeof(size_t)*number_of_timeseries);
        memset(p, 0, sizeof(size_t)*number_of_timeseries);
    }

    float *ts = (float*)malloc(sizeof(float) * length);

    int i, j, rep;
    for (i=1; i<=number_of_timeseries; i+=repetition)
    {
        generate(ts, length, r, normalize);

        for(rep=0; rep<repetition; rep++)
        {
            if (ists) {
                fwrite(ts, sizeof(float), length,data_file);
                if (hastimestamp) {
//                    long timestamp = 0;
                    time_t timestamp = time(nullptr);
                    std::cout << "time = " << timestamp << std::endl;
                    srand(i);
                    std::cout << "stime = " << rand() % timestamp << " " << i << std::endl;
                    long stimestamp = rand() % timestamp;
                    // 负数
                    fwrite(&stimestamp, sizeof(long), 1, data_file);
                }
            }
            if (issaxt) {
                saxt_from_ts(ts, saxts[i].asaxt);
                p[i] = i;

            }
//             for(j=0; j<length; j++) {
//                 printf ("%g ", ts[j]);
//             }
//             printf("\n");
        }
        if(i % (1000 * repetition) == 0) {
            fprintf(stderr,"\r\x1b[m>> Generating: \x1b[36m%2.2lf%%\x1b[0m",(float) ((float)i/(float)number_of_timeseries) * 100);
        }
    }
    fprintf(stderr, "\n");

    if (issaxt) {
        if (issort) {
            std::sort(saxts, saxts + number_of_timeseries + 1);
        }
        for (i=1; i<=number_of_timeseries; i+=repetition) {
            if (i%1000==0)
            saxt_print(saxts[i]);
            fwrite(p+i, sizeof(size_t), 1 , saxt_file);
            fwrite(saxts[i].asaxt, sizeof(saxt_type), Bit_cardinality , saxt_file);
        }
    }


    // Finalize random number generator
    if (ists) fclose (data_file);
    if (issaxt) fclose(saxt_file);
    gsl_rng_free (r);
}


//追加写入

int main(int argc, char **argv) {

    // Initialize variables
    int length = Ts_length;                 // Length of a single time series
    int number_of_timeseries = 3e6;   // Number of time series to generate

    // 倾斜度
    float skew_frequency = 0;           // The skew frequency
    int repetition = 1;             // How many times each time series is repeated
    char normalize = 1;             // Normalize or not.
    char * filename = "./output.bin";
    char * filename1 = "./saxt.bin";
    bool ists = 1;
    bool issaxt = 1;
    bool issort = 0;
    bool hastimestamp = 1;

//    // Parse command line arguments
//    parse_args(argc, argv, &length, &number_of_timeseries, &skew_frequency, &normalize,&filename);
//
//    fprintf(stderr,PRODUCT);
//
    if((1-skew_frequency) > 0)
        repetition = number_of_timeseries / (number_of_timeseries * (1-skew_frequency));
    else
        repetition = number_of_timeseries;
    fprintf(stderr, ">> Generating random time series...\n");
    fprintf(stderr, ">> Data Filename: %s\n", filename);
    generate_random_timeseries(length, number_of_timeseries, normalize, repetition,filename,issaxt,filename1,issort, ists, hastimestamp);
    fprintf(stderr, ">> Done.\n");
    return 0;
}