#include "helper.h"
#include "pfb.h"

// File containing helper functions for main

void printUsage(const char* progName) {
	(void) printf("Usage: %s [options] <data-file>\n", progName);
    (void) printf("    -h  --help                           ");
    (void) printf("Display this usage information\n");
    (void) printf("    -b  --nsub                           ");
    (void) printf("Number of sub-bands in the data\n");
    (void) printf("    -n  --nfft <value>                   ");
    (void) printf("Number of points in FFT\n");
	return;
}


int loadData(char* f, char* inputData) {
	int ret = EXIT_SUCCESS;
	int file =  0;

	int readSize = SAMPLES * DEF_NUM_CHANNELS * DEF_NUM_ELEMENTS * (2*sizeof(char));
	//inputData = (char*) malloc(readSize);
	if(NULL == inputData) {
		(void) fprintf(stderr, "ERROR: Memory allocation failed! %s.\n", strerror(errno));
		return EXIT_FAILURE;
	}

	file = open(f, O_RDONLY);
	if (file < EXIT_SUCCESS) {
		(void) fprintf(stderr, "ERROR: failed to open data file. %s\n", strerror(errno));
		return EXIT_FAILURE;
	}

	ret = read(file, inputData, readSize);
	if (ret < EXIT_SUCCESS) {
		(void) fprintf(stderr, "ERROR: failed to read data file. %s\n", strerror(errno));
		(void) close(file);
		return EXIT_FAILURE;
	}

	(void) close(file);
	return EXIT_SUCCESS;

}