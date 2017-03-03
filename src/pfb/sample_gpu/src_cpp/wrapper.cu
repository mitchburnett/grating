/*
* The wrapper.cu conatins the implementations of the functions declared in wrapper.h
*/

#include "wrapper.h"

// helper function defintions
void __checkCudaErrors(cudaError_t err, const char* const func, const char* file, const int line) {
	if(err != cudaSuccess) {
		std::cerr << "ERROR: file <" << file << ">" << ":" << line << " ";
		std::cerr << cudaGetErrorString(err) << " : " << func << std::endl;
		//cleanUp();
		//resetDevice();
		exit(0);
	}
}

void reduce(int* signal_d) {

	dim3 gridSize(1,1,1);
	dim3 blockSize(16,1,1);
	int smemsize = blockSize.x*sizeof(int);
	int n = 16;
	reduction<<<gridSize, blockSize, smemsize>>>(signal_d, n);
	checkCudaErrors(cudaGetLastError());

	return;
}