#include <stdlib.h>
#include <stdio.h>
#include <iostream>

//#include <cuda_runtime.h>

#include "wrapper.h"

using namespace std;

int main() {

	cout << "*** TESTING GPU WRAPPER IMPLEMENTATION ***" << endl;

	int* signal_h = (int *) malloc(10*sizeof(int));
	for(int i = 0; i < 10; i++) {
		signal_h[i] = i;
	}

	int* signal_d = NULL;
	checkCudaErrors(cudaMalloc((void **) &signal_d, 10*sizeof(int)));
	checkCudaErrors(cudaMemcpy(signal_d, signal_h, 10*sizeof(int), cudaMemcpyHostToDevice));

	reduce(signal_d);

	checkCudaErrors(cudaMemcpy(signal_h, signal_d, 1*sizeof(int), cudaMemcpyDeviceToHost));

	cout << "Reduction result: " << signal_h[0]  << endl;

	return 0;
}