#include "kernels.h"

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

__global__ void reduction(int* signal_d, int n) {

	int* smem = SharedMemory<int>();

	// load shared memory
	int tidx = threadIdx.x;
	int idx = blockIdx.x*blockDim.x + tidx;

	if(idx < n) {
		smem[tidx] = signal_d[idx];
	} else {
		smem[tidx] = 0;
	}

	__syncthreads();

	//perform reduction
	for(unsigned int s=blockDim.x/2; s > 0; s>>=1){
		if(tidx < s) {
			smem[tidx] += smem[tidx + s];
		}
		__syncthreads();
	}

	//get result and avg
	if(tidx == 0){
		signal_d[blockIdx.x] = smem[0];
	}

	return;
}