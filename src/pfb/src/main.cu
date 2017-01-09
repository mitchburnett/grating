#include "pfb.h"

char* g_inputData = NULL;
//char* g_inputData_d = NULL;
float2* g_outputData = NULL;

int loadData(char* f){
	int ret = EXIT_SUCCESS;
	int file =  0;

	int readSize = SAMPLES * DEF_NUM_CHANNELS * DEF_NUM_ELEMENTS * (2*sizeof(char));
	g_inputData = (char*) malloc(readSize);
	if(NULL == g_inputData) {
		(void) fprintf(stderr, "ERROR: Memory allocation failed! %s.\n", strerror(errno));
		return EXIT_FAILURE;
	}

	file = open(f, O_RDONLY);
	if (file < EXIT_SUCCESS) {
		(void) fprintf(stderr, "ERROR: failed to open data file. %s\n", strerror(errno));
		return EXIT_FAILURE;
	}

	ret = read(file, g_inputData, readSize);
	if (ret < EXIT_SUCCESS) {
		(void) fprintf(stderr, "ERROR: failed to read data file. %s\n", strerror(errno));
		(void) close(file);
		return EXIT_FAILURE;
	}

	(void) close(file);
	return EXIT_SUCCESS;

}

int main(int argc, char *argv[]) {

	int ret = EXIT_SUCCESS;

	/* valid short and long options */
	const char* const pcOptsShort = "hn:t:w:b:d:p";
	const struct option stOptsLong[] = {
		{ "help",		0, NULL,	'h' },   
		{ "nfft", 		1, NULL,	'n' },
		{ "taps",		1, NULL,	't' },
		{ "window",		1, NULL,	'w' },
		{ "nsub",		1, NULL,	'b' },
		{ "datatype",	1, NULL,	'd' },
		{ "plot",		0, NULL,	'p' },
		{ NULL,			0, NULL, 	0	}
	};

	const char* ProgName = argv[0];
	int argFlag = 0;

	/* parse input */
	int nextOpt = 0;

	// no arguments presented
	if(argc < optind) {
		(void) fprintf(stderr, "Missing required arguments\n");
		return EXIT_FAILURE;
	}

	// get data filename
	char filename[256] = {0};
	(void) strncpy(filename, argv[1], 256);
	filename[255] = '\0';

	// create coeff
	genCoeff(argc, argv);

	(void) fprintf(stdout, "Good Job!\n");
	return 0;
	// load data into memory
	ret = loadData(filename);
	if (ret == EXIT_FAILURE) {
		return EXIT_FAILURE;
	}

	// init cuda device
	int iCudaDevice = DEF_CUDA_DEVICE;
	ret = loadCoeff(iCudaDevice);

	// malloc data arrays
	//int inputSize = SAMPLES * DEF_NUM_CHANNELS * DEF_NUM_ELEMENTS * (2*sizeof(char));
	int outputSize = SAMPLES * PFB_CHANNELS * DEF_NUM_ELEMENTS * (2*sizeof(float)); // need to convince myself of this output data size.

	g_outputData = (float2*) malloc(outputSize);
	memset(g_outputData, 0, outputSize);

	// start pfb function
	int select = 0;
	ret = runPFB(g_inputData, g_outputData, select);
	if (ret == EXIT_FAILURE) {
		(void) fprintf(stderr, "ERROR: runPFB failed!\n");
		free(g_inputData);
		free(g_outputData);
		return EXIT_FAILURE;
	}

	// process return from pfb - write to file
	int file = 0;
	
	char outfile[256] = "output/outfile.dat\0";
	file = open(outfile,
					O_CREAT | O_TRUNC | O_WRONLY,
					S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	if(file < EXIT_SUCCESS) {
		(void) fprintf(stderr, "ERROR: writing outfile failed\n");
		free(g_inputData);
		free(g_outputData);
		return EXIT_FAILURE;
	}

	(void) write(file, g_outputData, outputSize);
	(void) close(file);

	free(g_inputData);
	free(g_outputData);

	return EXIT_SUCCESS;
}