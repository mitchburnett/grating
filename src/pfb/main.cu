#include "pfb.h"

int main(int argc, char *argv[]) {

	
	// setup data arrrays

	// load data into memory

	// init cuda device
	int iCudaDevice = DEF_CUDA_DEVICE;
	int iRet = loadCoeff(iCudaDevice);
	if(iRet != EXIT_SUCCESS) {
		(void) fprintf(stderr, "ERROR: Device init failed!\n");
		cleanUp();
		return EXIT_FAILURE;
	}

	// start pfb function

	// process return from pfb



	return EXIT_SUCCESS;
}