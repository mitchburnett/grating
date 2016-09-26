#include "pfb.h"

int main() {

	int iCudaDevice = DEF_CUDA_DEVICE;
	int iRet = loadCoeff(iCudaDevice);
	if(iRet != EXIT_SUCCESS) {
		(void) fprintf(stderr, "ERROR: Device init failed!\n");
		cleanUp();
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}