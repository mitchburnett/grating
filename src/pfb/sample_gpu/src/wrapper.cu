#include "wrapper.h"

void reduce(int* signal_d) {

	dim3 gridSize(1,1,1);
	dim3 blockSize(10,1,1);
	int smemsize = blockSize.x;
	int n = 10;
	reduction<<<gridSize, blockSize, smemsize>>>(signal_d, n);
	checkCudaErrors(cudaGetLastError());

	return;
}