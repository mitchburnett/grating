/* 
* Sample for building a c/c++ application using cuda. This sample shows
* how you would use a wrapper to seperate cuda specific syntax from c/c++
* syntax. 
*
*/

#include <stdlib.h>
#include <stdio.h>

#include "wrapper.h"

int main() {

	fprintf(stdout, "*** TESTING GPU WRAPPER IMPLEMENTATION ***\n");

	int* signal_h = (int *) malloc(16*sizeof(int));
	for(int i = 0; i < 16; i++) {
		signal_h[i] = i;
	}

	int* signal_d = NULL;
	checkCudaErrors(cudaMalloc((void **) &signal_d, 16*sizeof(int)));
	checkCudaErrors(cudaMemset((void *) signal_d, 0, 16*sizeof(int)));
	checkCudaErrors(cudaMemcpy(signal_d, signal_h, 16*sizeof(int), cudaMemcpyHostToDevice));

	reduce(signal_d);

	checkCudaErrors(cudaMemcpy(signal_h, signal_d, 16*sizeof(int), cudaMemcpyDeviceToHost));

	fprintf(stdout, "%d", signal_h[0]);

	free(signal_h);
	cudaFree(signal_d);

	return 0;
}