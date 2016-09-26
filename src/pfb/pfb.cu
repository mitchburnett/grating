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
//int g_iIsPFBOn = DEF_PFB_ON;

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
		int isPFB, int numSubBands, int nfft, int cudaDevice) {

	/*int iRet = EXIT_SUCCESS;
	int iSpecCount = 0;
	int NumAcc = DEF_ACC;
	*/
	//g_iIsPFBOn = isPFB;
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
	return 0;

}

// return true or false upon successful setup.
int loadCoeff(int iCudaDevice){

	int iRet = EXIT_SUCCESS;

	int iDevCount = 0;
	cudaDeviceProp stDevProp = {0};
	cufftResult iCUFFTRet = CUFFT_SUCCESS;
	size_t lTotCUDAMalloc = 0;

	int i = 0;

	//Register signal handlers?

	/********************************************/
	/* Look for eligable Cuda Device and select */
	/********************************************/
	(void) fprintf(stdout, "Querying CUDA devices.\n");

	(void) cudaGetDeviceCount(&iDevCount);
	if (0 == iDevCount) {
		(void) fprintf(stderr, "ERROR: No CUDA-capable device found!\n");
		return EXIT_FAILURE;
	}
	// Look for requested device (if applicable)
	if (iCudaDevice >= iDevCount) {
		(void) fprintf(stderr,
					   "ERROR: Requested device %d no found in present %d device list.\n",
					   iCudaDevice,
					   iDevCount);
		return EXIT_FAILURE;
	}
	// Query devices and setup selected device.
	for(i = 0; i < iDevCount; i++) {
		CUDASafeCallWithCleanUp(cudaGetDeviceProperties(&stDevProp, i));
		printf("\tDevice %d: %s, Compute Capability %d.%d, %d physical threads %s\n",
				i,
				stDevProp.name, stDevProp.major, stDevProp.minor,
				stDevProp.multiProcessorCount * stDevProp.maxThreadsPerMultiProcessor,
				(iCudaDevice == i) ? "<<SELECTED>>" : "");
	}
	CUDASafeCallWithCleanUp(cudaSetDevice(iCudaDevice));

	// Setup block and thread paramters
	CUDASafeCallWithCleanUp(cudaGetDeviceProperties(&stDevProp, 0));
	g_iMaxThreadsPerBlock = stDevProp.maxThreadsPerBlock;
	g_iMaxPhysThreads = stDevProp.multiProcessorCount * stDevProp.maxThreadsPerMultiProcessor;

	// Check if valid operation lengths. i.e. The input buffer is long enough (should this bee done here or elsewhere?)

	// Set malloc size - lTotCUDAMalloc is used only to calculate the total amount of memory not used for the allocation.
	lTotCUDAMalloc += g_iSizeRead; // size   data
	lTotCUDAMalloc += (g_iNumSubBands * g_iNFFT * sizeof(float(4))); // size of FFT input array This should be different since our data is unsigned char?
	lTotCUDAMalloc += (g_iNumSubBands * g_iNFFT * sizeof(float(4))); // size of FFT output array
	lTotCUDAMalloc += (g_iNumSubBands * g_iNFFT * sizeof(float)); 	// size of PFB Coefficients
	// Check CUDA device can handle the memory request
	if(lTotCUDAMalloc > stDevProp.totalGlobalMem) {
		(void) fprintf(stderr,
						"ERROR: Total memory requested on GPU is %g MB of %g possible MB.\n"
						"\tMemory break-down:\n"
						"\tInput data buffer:\t%g MB\n"
						"\tFFT in array:\t%g MB\n"
						"\tFFT out array:\t%g MB\n"
						"\tPFB Coefficients: %d KB\n",
						((float) lTotCUDAMalloc) / (1024*1024),
						((float) stDevProp.totalGlobalMem) / (1024*1024),
						((float) g_iSizeRead) / (1024 * 1024),
						((float) g_iNumSubBands * g_iNFFT * sizeof(float4)) / (1024 * 1024),
						((float) g_iNumSubBands * g_iNFFT * sizeof(float4)) / (1024 * 1024)),
						((float) g_iNumSubBands * g_iNFFT * sizeof(float));
		return EXIT_FAILURE;
	}
	
	// print memory usage report.
	(void) fprintf(stderr,
					"INFO: Total memory requested on GPU is %g MB of %g possible MB.\n"
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

	/*************************/
	/* Load PFB coefficients */
	/*************************/
	(void) fprintf(stdout, "Setting up PFB filter coefficients...\n");
	g_iNTaps = NUM_TAPS; // set the number of taps. Change this to where it happens earlier to be more dynamic.
	int sizePFB = g_iNumSubBands * g_iNTaps * g_iNFFT * sizeof(float);

	// Allocate memory for PFB coefficients to be read in
	g_pfPFBCoeff = (float *) malloc(sizePFB); // allocate the memory needed for the size of one pfb pass through
	if(NULL == g_pfPFBCoeff) {
		(void) fprintf(stderr, "ERROR: Memory allocation for the PFB coefficients failed. %s\n",
								strerror(errno));
		return EXIT_FAILURE;
	}

	// Read filter coefficients from file
	(void) fprintf(stdout, "Reading in coefficients...\n");
	(void) sprintf(g_acFileCoeff,
				   "%s_%s_%d_%d_%d%s",
				   FILE_COEFF_PREFIX,
				   FILE_COEFF_DATATYPE,
				   g_iNTaps,
				   g_iNFFT,
				   g_iNumSubBands,
				   FILE_COEFF_SUFFIX);

	g_iFileCoeff = open(g_acFileCoeff, O_RDONLY);
	if(g_iFileCoeff < EXIT_SUCCESS) {
		(void) fprintf(stderr, "ERROR: Failed to open coefficient file %s. %s\n",
					  			g_acFileCoeff,
					  			strerror(errno));
		return EXIT_FAILURE;
	}

	iRet = read(g_iFileCoeff, g_pfPFBCoeff, sizePFB);
	if(iRet != sizePFB) {
		(void) fprintf(stderr, "ERROR: Failed reading filter coefficients. %s\n", strerror(errno));
		return EXIT_FAILURE;
	}
	(void) close(g_iFileCoeff);

	/********************************************/
	/* Allocate memory and setup on CUDA device */
	/********************************************/
	(void) fprintf(stdout, "Setting up CUDA device.\n");

	// allocate memory for pfb coefficients on GPU
	(void) fprintf(stdout, "\tAllocating memory for PFB...\n");
	CUDASafeCallWithCleanUp(cudaMalloc((void **) &g_pfPFBCoeff_d, sizePFB));

	// copy coeff to device
	(void) fprintf(stdout, "\tCopying filter coefficients...\n");
	CUDASafeCallWithCleanUp(cudaMemcpy(g_pfPFBCoeff_d, g_pfPFBCoeff, sizePFB, cudaMemcpyHostToDevice));

	// allocate memory for FFT in and out arrays
	(void) fprintf(stdout, "\tAllocate memory for FFT arrays...\n");
	int sizeDataBlock = g_iNumSubBands * g_iNFFT * sizeof(float4);
	CUDASafeCallWithCleanUp(cudaMalloc((void **) &g_pf4FFTIn_d, sizeDataBlock));
	CUDASafeCallWithCleanUp(cudaMalloc((void **) &g_pf4FFTOut_d, sizeDataBlock));

	// set kernel parameters
	(void) fprintf(stdout, "\tSetting kernel parameters...\n");
	if(g_iNFFT < g_iMaxThreadsPerBlock) {
		g_dimBPFB.x   = g_iNFFT;
		g_dimBCopy.x = g_iNFFT;
	} else {
		g_dimBPFB.x   = g_iMaxThreadsPerBlock;
		g_dimBCopy.x = g_iMaxThreadsPerBlock;
	}
	g_dimGPFB.x  = (g_iNumSubBands * g_iNFFT) / g_dimGPFB.x;
	g_dimGCopy.x = (g_iNumSubBands * g_iNFFT) / g_dimBCopy.x;

	// create a CUFFT plan
	(void) fprintf(stdout, "\tCreating cuFFT plan...\n");
	iCUFFTRet = cufftPlanMany(&g_stPlan,
							  FFTPLAN_RANK,
							  &g_iNFFT,
							  &g_iNFFT,
							  FFTPLAN_ISTRIDE,
							  FFTPLAN_IDIST,
							  &g_iNFFT,
							  FFTPLAN_OSTRIDE,
							  FFTPLAN_ODIST,
							  CUFFT_C2C,
							  FFTPLAN_BATCH);
	if(iCUFFTRet != CUFFT_SUCCESS) {
		(void) fprintf(stderr, "ERROR: Plan creation failed!\n");
		return EXIT_FAILURE;
	}

	fprintf(stdout, "\nDevice for PFB successful initialized!\n");
	return EXIT_SUCCESS;

}

void __CUDASafeCallWithCleanUp(cudaError_t iRet,
                               const char* pcFile,
                               const int iLine,
                               void (*pcleanUp)(void))
{
    if (iRet != cudaSuccess)
    {
        (void) fprintf(stderr,
                       "ERROR: File <%s>, Line %d: %s\n",
                       pcFile,
                       iLine,
                       cudaGetErrorString(iRet));
        /* free resources */
        (*pcleanUp)();
        exit(EXIT_FAILURE);
    }

    return;
}

void cleanUp() {
/* free resources */
    if (g_pc4InBuf != NULL) {
        free(g_pc4InBuf);
        g_pc4InBuf = NULL;
    }
    if (g_pc4Data_d != NULL) {
        (void) cudaFree(g_pc4Data_d);
        g_pc4Data_d = NULL;
    }
    if (g_pf4FFTIn_d != NULL) {
        (void) cudaFree(g_pf4FFTIn_d);
        g_pf4FFTIn_d = NULL;
    }
    if (g_pf4FFTOut_d != NULL) {
        (void) cudaFree(g_pf4FFTOut_d);
        g_pf4FFTOut_d = NULL;
    }

    free(g_pfPFBCoeff);
    (void) cudaFree(g_pfPFBCoeff_d);

    /* destroy plan */
    /* TODO: check for plan */
    (void) cufftDestroy(g_stPlan);

    return;
}















