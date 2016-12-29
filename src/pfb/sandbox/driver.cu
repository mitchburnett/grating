//#include "driver.cu"

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cufft.h>

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

#define NUM_EL 		 64
#define CHANNELS 	 25
#define PFB_CHANNELS 5
#define SAMPLES		 4000

#define DEF_CUDA_DEVICE 0

#define checkCudaErrors(err) __checkCudaErrors(err, __FILE__, __LINE__)


char* g_inputData = NULL;
char2* g_outputData = NULL;
char* g_inputData_d = NULL;
char2* g_outputData_d = NULL;

int loadData(char* f){
	int ret = EXIT_SUCCESS;
	int file =  0;

	int readSize = NUM_EL * CHANNELS * SAMPLES * (2*sizeof(char));
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

void __checkCudaErrors(cudaError_t err, const char* file, const int line) {
	if (err != cudaSuccess) {
		(void) fprintf(stderr, "ERROR: file <%s>, Line %d: %s\n",
						file,
						line,
						cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}  

	return;
}

int init(){

	int cudaDevice = DEF_CUDA_DEVICE;
	checkCudaErrors(cudaSetDevice(cudaDevice));

	int inputSize  = NUM_EL * CHANNELS * SAMPLES * (2*sizeof(char));
	int outputSize = SAMPLES * PFB_CHANNELS * NUM_EL * (2*sizeof(char));

	// allocate memory for input and output data on the device.
	checkCudaErrors(cudaMalloc((void **) &g_inputData_d, inputSize));
	checkCudaErrors(cudaMemset((void *) g_inputData_d, 0, inputSize));
	checkCudaErrors(cudaMalloc((void **) &g_outputData_d, outputSize));
	checkCudaErrors(cudaMemset((void *) g_outputData_d, 0, outputSize));

	// copy data to the device.
	checkCudaErrors(cudaMemcpy(g_inputData_d, g_inputData, inputSize, cudaMemcpyHostToDevice));

	return EXIT_SUCCESS;
}

__global__ void map(char* dataIn,
			   		char2* dataOut,
			   		int channelSelect) {

	// select the channel range
	int channelMin = PFB_CHANNELS*channelSelect;
	
	int absIdx = 2 * blockDim.y*(blockIdx.x*CHANNELS + (channelMin+blockIdx.y)) + 2 * threadIdx.y; // times 2 because we are mapping a sequence of values to char2 array.
	int mapIdx = blockDim.y*(blockIdx.x*gridDim.y + blockIdx.y) + threadIdx.y;

	dataOut[mapIdx].x = dataIn[absIdx];
	dataOut[mapIdx].y = dataIn[absIdx+1];
	return;
}

int main(int argc, char *argv[]) {

	int ret = EXIT_SUCCESS;
	if(argc < 2) {
		(void) fprintf(stderr, "ERROR: Data filename not specified.\n");
		return EXIT_FAILURE;
	}

	char filename[256] = {0};
	(void) strncpy(filename, argv[1], 256);
	filename[255] = '\0';

	ret = loadData(filename);
	if (ret == EXIT_FAILURE) {
		return EXIT_FAILURE;
	}

	ret = init();

	// run map
	int select = 0;
	dim3 gridSize(SAMPLES,PFB_CHANNELS,1);
	dim3 blockSize(1, NUM_EL, 1);
	map<<<gridSize, blockSize>>>(g_inputData_d, g_outputData_d, select);
	checkCudaErrors(cudaGetLastError());	



	int outputSize = SAMPLES * PFB_CHANNELS * NUM_EL * (2*sizeof(char));
	g_outputData = (char2*) malloc(outputSize);
	checkCudaErrors(cudaMemcpy(g_outputData, g_outputData_d, outputSize, cudaMemcpyDeviceToHost));

	//output the true data as a check.
	/*int file = 0;
	char outfileFull[256] = "outfileFull.dat\0";
	file = open(outfile,
					O_CREAT | O_TRUNC | O_WRONLY,
					S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	if(file < EXIT_SUCCESS) {
		(void) fprintf(stderr, "ERROR: writing outfile failed\n");
		return EXIT_FAILURE;
	}

	(void) write(file, g_inputData, SAMPLES*CHANNELS*NUM_EL*2*sizeof(char));
	(void) close(file); */

	
	// output the mapped data.
	int file = 0;
	char outfile[256] = "outfile.dat\0";
	file = open(outfile,
					O_CREAT | O_TRUNC | O_WRONLY,
					S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	if(file < EXIT_SUCCESS) {
		(void) fprintf(stderr, "ERROR: writing outfile failed\n");
		return EXIT_FAILURE;
	}

	(void) write(file, g_outputData, outputSize);
	(void) close(file);

	return EXIT_SUCCESS;

}








