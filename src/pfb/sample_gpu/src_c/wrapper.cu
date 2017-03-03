/*
* The wrapper.cu conatins the implementations of the functions declared in wrapper.h
*/
#ifdef __cplusplus
extern "C" {
#include "wrapper.h"
}
#endif
// helper function defintions
void __checkCudaErrors(cudaError_t err, const char* const func, const char* file, const int line) {
	if(err != cudaSuccess) {
		fprintf(stderr, "ERROR: file <%s> : %d\n", file, line );
		fprintf(stderr, "%s : %s\n",cudaGetErrorString(err), func);
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