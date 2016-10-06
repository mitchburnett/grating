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
#define SCALE_FACTOR  127
#define F_S		      256   // MHz
#define N			  4000  // Time samples
#define CHANNELS	  25    // Freq Channels
#define NUM_EL		  64    // Antenna Elements

int main(int argc, char *argv[]) {

	if(argc < 2) {
		(void) fprintf(stderr, "ERROR: Data filename not specifie.\n");
		return EXIT_FAILURE;
	}

	int iFile = 0;
	char acDataFilename[LEN_GENSTRING] = {0};
	(void) strncpy(acDataFilename, argv[1], LEN_GENSTRING);
	acDataFilename[LEN_GENSTRING-1] = '\0'; //NUll terminator at end of filename string.

	iFile = open(acDataFilename,
					O_CREAT | O_TRUNC | O_WRONLY,
					S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	if(EXIT_FAILURE == iFile) {
		(void) fprintf(stderr, "ERROR: Failed to open output file. %s\n", strerror(errno));
		return EXIT_FAILURE;
	}

	int i = 0;

	// generate freq array of
	int channelBandgap = 5;		// MHz jumps
	float freq[CHANNELS] = {};
	for(i = 1; i <= CHANNELS; i++) {
		freq[i] = channelBandgap * i;
		fprintf(stdout, "Freq Channel: %f\n", freq[i]);
	}

	int n = 0;
	int f = 0;
	int e = 0;
 
	signed char cDataReX = 0;
	signed char cDataImY = 0;
	signed char toWrite[NUM_EL*CHANNELS*(2*sizeof(char))] = {};
	int size = NUM_EL*CHANNELS*(2*sizeof(char));
	for(n = 0; n < N; n++) {
		for(f = 0; f < CHANNELS; f++) {

			//use the same sample for all elements
			cDataReX = SCALE_FACTOR * (0.1 * cos(2*M_PI * freq[f] * n / F_S));
			cDataImY = SCALE_FACTOR * (0.1 * sin(2*M_PI * freq[f] * n / F_S));
			for(e = 0; e < 2*NUM_EL; e++) {
				
				int idx = e + f * NUM_EL;
				if( !(e%2) ) {
				//create interleaved samples for real and Im 
				toWrite[idx] = cDataReX;
				} else {
				toWrite[idx] = cDataImY;
				}
			}
		}
		(void) write(iFile, toWrite, NUM_EL*CHANNELS*(2*sizeof(char)));
	}

	(void) close(iFile);

	return EXIT_SUCCESS;
}