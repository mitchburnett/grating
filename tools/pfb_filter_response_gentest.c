//
//  pfb_filter_response_gentest.c
//  pfb_filter_response_gentest
//
//  Created by Mitchell Burnett on 9/6/16.
//  Copyright © 2016 Mitchell Burnett. All rights reserved.
//

#include <stdio.h>

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <float.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <assert.h>


#define LEN_GENSTRING 256
#define DEF_NUM_SAMPS 256       /* Default number of samples */
#define F_S           256.0     /* sampling frequency in MHz */

#define SCALE_FACTOR  127       /* scale factor to keep values in char range */

int main(int argc, char *argv[]) {
    
    signed char cDataReX[DEF_NUM_SAMPS] = {};
    signed char cDataImX[DEF_NUM_SAMPS] = {};
    signed char cDataReY[DEF_NUM_SAMPS] = {};
    signed char cDataImY[DEF_NUM_SAMPS] = {};
    
    int i = 0;
    int j = 0;
    
    int iSampleLength = DEF_NUM_SAMPS;
    int iFile = 0;
    char acFileData[LEN_GENSTRING] = {0};
    
    float afFreqX[256] = {};
    float afFreqY[256] = {};
    
    const char *pcProgName = NULL;
    
    //Declare command line options.
    
    //Parse command line inputs.
    if(argc < 2) {
        (void) fprintf(stderr, "ERROR: Data file not specified!\n");
        return EXIT_FAILURE;
    }
    
    (void) strncpy(acFileData, argv[1], LEN_GENSTRING);
    acFileData[LEN_GENSTRING-1] = '\0'; //Null terminator at end of file string.
    
    iFile = open(acFileData,
                 O_CREAT | O_TRUNC | O_WRONLY,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if(EXIT_FAILURE == iFile) {
        (void) fprintf(stderr, "ERROR: Opening output file failed! %s.\n", strerror(errno));
        return EXIT_FAILURE;
    }

    //Generate freq array.
    for(i = 1; i <= 256; i++) {
        afFreqX[i-1] = i*1.0;
        afFreqY[i-1] = afFreqX[i-1];
    }
    
    for(i = 0; i < 256; i++) {
        for(j = 0; j < 256; j++) {
        cDataReX[j] = SCALE_FACTOR*(0.1*cos(2*M_PI *afFreqX[i] * j /F_S));
        cDataImX[j] = SCALE_FACTOR*(0.1*sin(2*M_PI *afFreqX[i] * j /F_S));
        }
        (void) write(iFile, &cDataReX, 256*sizeof(signed char));
        (void) write(iFile, &cDataImX, 256*sizeof(signed char));
    }
    
    for(i = 0; i < 256; i++) {
        for(j = 0; j < 256; j++) {
            cDataReY[j] = SCALE_FACTOR*(0.1*cos(2*M_PI *afFreqY[i] * j /F_S));
            cDataImY[j] = SCALE_FACTOR*(0.1*sin(2*M_PI *afFreqY[i] * j /F_S));
        }
        (void) write(iFile, &cDataReY, 256*sizeof(signed char));
        (void) write(iFile, &cDataImY, 256*sizeof(signed char));
    }
    
    (void) close(iFile);
    
    return EXIT_SUCCESS;
}











