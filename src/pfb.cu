#include "pfb.h"

/*

dim3 g_dimBAccum(1, 1, 1);
dim3 g_dimGAccum(1, 1);

float4* g_pf4SumStokes = NULL;
float4* g_pf4SumStokes_d = NULL;

char g_acFileData[256] = {0}; // File data to load and process. If this is a function data is an input.

int g_iNumSubBands = DEF_NUM_SUBBANDS;
*/

int g_IsDataReadDone = FALSE;
int g_IsProcDone = FALSE;
int g_iIsPFBOn = DEF_PFB_ON;

int g_iSizeFile = 0;
int g_iReadCount = 0;
int g_iSizeRead = DEF_SIZE_READ;
int g_iFileCoeff = 0;
char g_acFileCoeff[256] = {0};

int g_iMaxThreadsPerBlock = 0;
int g_iMaxPhysThreads;
dim3 g_dimBPFB(1, 1, 1);
dim3 g_dimGPFB(1, 1);
dim3 g_dimBCopy(1, 1, 1);
dim3 g_dimGCopy(1, 1);
cufftHandle g_stPlan = {0};

char4* g_pc4InBuf = NULL;
char4* g_pc4InBufRead = NULL;
char4* g_pc4Data_d = NULL;
char4* g_pc4DataRead_d = NULL;
float4* g_pf4FFTIn_d = NULL;
float4* g_pf4FFTOut_d = NULL;

int g_iNFFT = DEF_LEN_SPEC;
int g_iNTaps = NUM_TAPS;
int g_iNumSubBands = DEF_NUM_SUBBANDS;
float *g_pfPFBCoeff = NULL;
float *g_pfPFBCoeff_d = NULL;

// The main will potentially be a different function if this is part of a library?
// inputs: numSubbands, nfft, isPFBOn, iCudaDevice
int pfb(unsigned char* inputData_h,
		unsigned char* outputData_h,
		int, isPFB, int numSubBands, int nfft, int cudaDevice) {

	/*int iRet = EXIT_SUCCESS;
	int iSpecCount = 0;
	int NumAcc = DEF_ACC;
	*/
	g_iIsPFBOn = isPFB;
	g_iNFFT = nfft;
	g_iNumSubBands = numSubBands;
	int iCUDADevice = cudaDevice;

	cudaError_t iCUDARet = cudaSuccess;

	int iProcData = 0;
	long int lProcDataAll = 0;

	// Time vars without deep benchmarking
	struct timeval stStart = {0};
	struct timeval stStop = {0};
	float fTimeTaken = 0.0;
	float fTotThroughput = 0.0;

	/* Init */

}

// return true or false upon successful setup.
int loadCoeff(int iCudaDevice){

	int iRet = EXIT_SUCCESS;

	int iDevCount = 0;
	cudaDeviceProp stDevProp = {0};
	cufftResult iCUFFTRet = CUFFT_SUCCES;
	size_t lTotCUDAMalloc = 0;

	int i = 0;

	//Register signal handlers?

	// Look for Cuda Devices
	(void) cudaGetDeviceCount(&iDevCount);
	if (0 == iDevCount) {
		(void) fprintf(stder, "ERROR: No CUDA-capable device found!\n");
		return EXIT_FAILURE;
	}
	// Loof for requested device (if applicable)
	if (iCUDADevice >= iDevCount) {
		(void) fprintf(stderr,
					   "ERROR: Requested device %d no found in present %d device list.\n",
					   iCUDADevice,
					   iDevCount);
		return EXIT_FAILURE;
	}
	// Query devices and setup selected device.
	for(i = 0; i < iDevCount; i++) {
		CUDASafeCallWithCleanUp(cudaGetDeviceProperties(&stDevProp, i));
		printf("Device %d: %s, Compute Capability %d.%d, %d physical threads %s\n",
				i,
				stDevProp.name, stDevProp.major, stDevProp.minor,
				stDevProp.multiProcessorCount * stDevProp.maxThreadsPerMultiProcessor,
				(iCUDADevice == i) ? "selected" : "");
	}
	CUDASafeCallWithCleanUp(cudaSetDevice(iCudaDevice));

	// Setup block and thread paramters
	CUDASafeCallWithCleanUp(cudaGetDeviceProperties(&stDevProp, 0));
	g_iMaxThreadsPerBlock = stDevProp.maxThreadsPerBlock;
	g_iMaxPhysThreads = stDevProp.multiProcessorCount * stDevProp.maxThreadsPerMultiProcessor;

	// Check if valid operation lengths. i.e. The input buffer is long enough (should this bee done here or elsewhere?)

	// Set malloc size - lTotCUDAMalloc is used only to calculate the total amount of memory not used for the allocation.
	lTotCUDAMalloc += g_iSizeRead; // size   data
	lTotCUDAMalloc += (g_iNumSubBands * g_iNFFT * sizeof(float(4))) // size of FFT input array This should be different since our data is unsigned char?
	lTotCUDAMalloc += (g_iNumSubBands * g_iNFFT * sizeof(float(4))) // size of FFT output array
	lTotCUDAMalloc += (g_iNumSubBands * g_iNFFT * sizeof(float)) 	// size of PFB Coefficients
	// Check CUDA device can handle the memory request
	if(lTotCUDAMalloc > stDevProp.totalGlobalMem) {
		(void) fprintf(stderr,
						"ERROR: Total memory requested on GPU is %g MB of %g possible MB.\n"
						"Memory break-down:\n"
						"\tInput data buffer:\t%g MB\n"
						"\tFFT in array:\t%g MB"
						"\tFFT out array:\t%g MB"
						"\tPFB Coefficients: %d KB\n",
						((float) lTotCUDAMalloc) / (1024*1024),
						((float) stDevProp.totalGlobalMem) / (1024*1024),
						((float) g_iSizeRead) / (1024 * 1024),
						((float) g_iNumSubBands * g_iNFFT * sizeof(float4)) / (1024 * 1024),
						((float) g_iNumSubBands * g_iNFFT * sizeof(float4)) / (1024 * 1024)),
						((float) g_iNumSubBands * g_iNFFT * sizeof(float));
	}
#ifdef DEBUG
	(void) fprintf(stderr,
					"ERROR: Total memory requested on GPU is %g MB of %g possible MB.\n"
					"Memory break-down:\n"
					"\tInput data buffer:\t%g MB\n"
					"\tFFT in array:\t%g MB"
					"\tFFT out array:\t%g MB"
					"\tPFB Coefficients: %d KB\n",
					((float) lTotCUDAMalloc) / (1024*1024),
					((float) stDevProp.totalGlobalMem) / (1024*1024),
					((float) g_iSizeRead) / (1024 * 1024),
					((float) g_iNumSubBands * g_iNFFT * sizeof(float4)) / (1024 * 1024),
					((float) g_iNumSubBands * g_iNFFT * sizeof(float4)) / (1024 * 1024)),
					((float) g_iNumSubBands * g_iNFFT * sizeof(float));
#endif



}	


















