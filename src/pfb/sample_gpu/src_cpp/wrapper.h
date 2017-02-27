/*
* The wrapper serves as the handshake between cuda and c/c++. This includes the declarations
* of the functions that are used to interface and run the cuda kernels.
*
*/

#include <stdio.h>
extern "C" {
	#include "kernels.h"
}

#define checkCudaErrors(err) __checkCudaErrors(err, #err, __FILE__, __LINE__)
void __checkCudaErrors(cudaError_t err, const char* const func, const char* file, const int line);

void reduce(int* signal_d);