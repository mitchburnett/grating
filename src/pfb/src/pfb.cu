#include "pfb.h"

/*

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

char2* g_pc2InBuf = NULL;
char2* g_pc2InBufRead = NULL;

char2* g_pc2Data_d = NULL;
char2* g_pc2DataRead_d = NULL;

float2* g_pf2FFTIn_d = NULL;
float2* g_pf2FFTOut_d = NULL;

int g_iNFFT = DEF_LEN_SPEC;
int g_iNTaps = NUM_TAPS;
//int g_iNumSubBands = DEF_NUM_CHANNELS * DEF_NUM_ELEMENTS;
int g_iNumSubBands = PFB_CHANNELS * DEF_NUM_ELEMENTS;

float *g_pfPFBCoeff = NULL;
float *g_pfPFBCoeff_d = NULL;

char* g_pcInputData_d = NULL;

int runPFB(char* inputData_h,
		   float2* outputData_h,
		   params pfbParams) {

	int channelSelect = pfbParams.select;

	//process variables
	int iRet = EXIT_SUCCESS;
	int countPFB = 0; // count number of times pfb fires.
	int countCpyFFT = 0;
	int countFFT = 0; // count number of FFT's computed.
	long lProcData = 0; // count how much data processed
	long ltotData = SAMPLES * PFB_CHANNELS * DEF_NUM_ELEMENTS; // total amount of data to proc

	//malloc and copy data to device
	int fullSize = SAMPLES * DEF_NUM_CHANNELS * DEF_NUM_ELEMENTS * (2*sizeof(char));
	int mapSize = SAMPLES * PFB_CHANNELS * DEF_NUM_ELEMENTS * (2*sizeof(char));
	CUDASafeCallWithCleanUp(cudaMalloc((void **) &g_pcInputData_d, fullSize));
	CUDASafeCallWithCleanUp(cudaMemset((void *)   g_pcInputData_d, 0, fullSize));
	CUDASafeCallWithCleanUp(cudaMalloc((void **) &g_pc2Data_d, mapSize));
	CUDASafeCallWithCleanUp(cudaMemset((void *)   g_pc2Data_d, 0, mapSize));

	CUDASafeCallWithCleanUp(cudaMemcpy(g_pcInputData_d, inputData_h, fullSize, cudaMemcpyHostToDevice));

	// extract channel data from full data stream and load into buffer.
	dim3 mapGSize(SAMPLES, PFB_CHANNELS, 1);
	dim3 mapBSize(1, DEF_NUM_ELEMENTS, 1);
	map<<<mapGSize, mapBSize>>>(g_pcInputData_d, g_pc2Data_d, channelSelect);
	CUDASafeCallWithCleanUp(cudaGetLastError());
	CUDASafeCallWithCleanUp(cudaThreadSynchronize());

	//Arrays for debugging the output structure.
	float2* fftIn = NULL;
	fftIn = (float2*) malloc(SAMPLES*PFB_CHANNELS*DEF_NUM_ELEMENTS*sizeof(float2));
	memset(fftIn, 0, SAMPLES*PFB_CHANNELS*DEF_NUM_ELEMENTS*sizeof(float2));
	// Begin PFB

	// p_pc2Data_d contains all the data. DataRead will update with each pass through the PFB.
	g_pc2DataRead_d = g_pc2Data_d;
	int pfb_on = 1;
	while(!g_IsProcDone){

		if(pfb_on) {
			//PFB
			PFB_kernel<<<g_dimGPFB, g_dimBPFB>>>(g_pc2DataRead_d, g_pf2FFTIn_d, g_pfPFBCoeff_d);
			CUDASafeCallWithCleanUp(cudaGetLastError());
			CUDASafeCallWithCleanUp(cudaThreadSynchronize());

			//update data read pointer
			g_pc2DataRead_d += g_iNumSubBands * g_iNFFT;
			++countPFB;
		} else {
			CopyDataForFFT<<<g_dimGPFB, g_dimBPFB>>>(g_pc2DataRead_d, g_pf2FFTIn_d);

			g_pc2DataRead_d += g_iNumSubBands * g_iNFFT;
			++countCpyFFT;
		}

		//copy pre fft in data to compare output with straigh fft.
		CUDASafeCallWithCleanUp(cudaMemcpy(fftIn, g_pf2FFTIn_d, g_iNumSubBands*g_iNFFT*sizeof(float2), cudaMemcpyDeviceToHost));
		fftIn += g_iNumSubBands * g_iNFFT;

		//FFT
		iRet = doFFT();
		if(iRet != EXIT_SUCCESS) {
			(void) fprintf(stderr, "ERROR: FFT failed\n");
			cleanUp();
			return EXIT_FAILURE;
		}
		CUDASafeCallWithCleanUp(cudaGetLastError());
		++countFFT;

		// copy data back to host.
		int outDataSize = g_iNumSubBands * g_iNFFT * (sizeof(cufftComplex));
		CUDASafeCallWithCleanUp(cudaMemcpy(outputData_h, g_pf2FFTOut_d, outDataSize, cudaMemcpyDeviceToHost));

		//update output data pointer.
		outputData_h += g_iNumSubBands * g_iNFFT;

		//update proc data
		lProcData += g_iNumSubBands * g_iNFFT;
		(void) fprintf(stdout, "Counters--PFB:%d FFT:%d\n",countPFB, countFFT);
		(void) fprintf(stdout, "Data process by the numbers:\n Processed:%ld\n To Process:%ld\n\n",lProcData, ltotData);
		if(lProcData > ltotData - NUM_TAPS*g_iNumSubBands*g_iNFFT){
			g_IsProcDone = TRUE;
		}

	}
	// reset pointer
	fftIn -= g_iNumSubBands*g_iNFFT*117;
	int file = 0;
	
	char outfile[256] = "output/prefft_out.dat\0";
	file = open(outfile,
					O_CREAT | O_TRUNC | O_WRONLY,
					S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	if(file < EXIT_SUCCESS) {
		(void) fprintf(stderr, "ERROR: writing outfile failed\n");
		cleanUp();
		resetDevice();
		return EXIT_FAILURE;
	}
	(void) write(file, fftIn, SAMPLES*PFB_CHANNELS*DEF_NUM_ELEMENTS*sizeof(float2));

	cleanUp();
	iRet = resetDevice();

	return iRet;

}

// make a call to execute a ptyhon program.
void genCoeff(int argc, char* argv[], params pfbParams) {

	FILE* file;
	char fname[256] = {"../../../scripts/grating_gencoeff.py"};

	int argCount = 11;
	char* arguments[32] = {}; // come back and create a dynamic structure, i.e definetly do not need 32, always 10 or 11.
	int i = 0;
	for(i = 0; i < 32; i++) {
		arguments[i] = (char*) malloc(256*sizeof(char*));
	}

	arguments[0] = argv[0];

	arguments[1] = (char*) "-n\0"; // (char*) acknowledges that I am assigning a const literal to a mutable and removes compile warnings for now.
	sprintf(arguments[2], "%d", pfbParams.nfft);

	arguments[3] = (char*) "-t\0";
	sprintf(arguments[4], "%d", pfbParams.taps);

	arguments[5] = (char*) "-b\0";
	sprintf(arguments[6], "%d", pfbParams.subbands);

	arguments[7] = (char*) "-w\0";
	sprintf(arguments[8], "%s", pfbParams.window);

	arguments[9] = (char*) "-d\0";
	sprintf(arguments[10], "%s", pfbParams.dataType);

	if(pfbParams.plot) {
		arguments[11] = (char*) "-p\0";
		argCount++;
	}

	for(i = 0; i < argCount; i++){
		fprintf(stdout, " %s", arguments[i]);
	}
	fprintf(stdout, "\n");

	// initalize and run python script
	Py_SetProgramName(argv[0]);
	Py_Initialize();
	PySys_SetArgv(argCount, arguments);
	file = fopen(fname, "r");
	PyRun_SimpleFile(file, fname);
	Py_Finalize();

	return;
}

// return true or false upon successful setup.
int initPFB(int iCudaDevice, params pfbParams){

	int iRet = EXIT_SUCCESS;

	// set pfb params from input parameters.
	g_iNFFT = pfbParams.nfft;
	g_iNTaps = pfbParams.taps;

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
	lTotCUDAMalloc += (g_iNumSubBands * g_iNFFT * sizeof(float(2))); // size of FFT input array This should be different since our data is unsigned char?
	lTotCUDAMalloc += (g_iNumSubBands * g_iNFFT * sizeof(float(2))); // size of FFT output array
	lTotCUDAMalloc += (g_iNumSubBands * g_iNFFT * sizeof(float)); 	// size of PFB Coefficients
	// Check CUDA device can handle the memory request
	if(lTotCUDAMalloc > stDevProp.totalGlobalMem) {
		(void) fprintf(stderr,
						"ERROR: Total memory requested on GPU is %g MB of %g possible MB.\n"
						"\t**** Memory breakdown *****\n"
						"\tInput data buffer:\t%g MB\n"
						"\tFFT in array:\t%g MB\n"
						"\tFFT out array:\t%g MB\n"
						"\tPFB Coefficients: %f KB\n",
						((float) lTotCUDAMalloc) / (1024*1024),
						((float) stDevProp.totalGlobalMem) / (1024*1024),
						((float) g_iSizeRead) / (1024 * 1024),
						((float) g_iNumSubBands * g_iNFFT * sizeof(float2)) / (1024 * 1024),
						((float) g_iNumSubBands * g_iNFFT * sizeof(float2)) / (1024 * 1024),
						((float) g_iNumSubBands * g_iNFFT * sizeof(float)));
		return EXIT_FAILURE;
	}
	
	// print memory usage report.
	(void) fprintf(stdout,
					"INFO: Total memory requested on GPU is %g MB of %g possible MB.\n"
					"\t**** Memory breakdown ****\n"
					"\tInput data buffer:\t%g MB\n"
					"\tFFT in array:\t%g MB\n"
					"\tFFT out array:\t%g MB\n"
					"\tPFB Coefficients: %f KB\n",
					((float) lTotCUDAMalloc) / (1024*1024),
					((float) stDevProp.totalGlobalMem) / (1024*1024),
					((float) g_iSizeRead) / (1024 * 1024),
					((float) g_iNumSubBands * g_iNFFT * sizeof(float2)) / (1024 * 1024),
					((float) g_iNumSubBands * g_iNFFT * sizeof(float2)) / (1024 * 1024),
					((float) g_iNumSubBands * g_iNFFT * sizeof(float)));

	/*************************/
	/* Load PFB coefficients */
	/*************************/
	(void) fprintf(stdout, "\nSetting up PFB filter coefficients...\n");
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
	(void) fprintf(stdout, "\tReading in coefficients...\n");
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
	(void) fprintf(stdout, "\nSetting up CUDA device.\n");

	// allocate memory for pfb coefficients on GPU
	(void) fprintf(stdout, "\tAllocating memory for PFB...\n");
	CUDASafeCallWithCleanUp(cudaMalloc((void **) &g_pfPFBCoeff_d, sizePFB));

	// copy coeff to device
	(void) fprintf(stdout, "\tCopying filter coefficients...\n");
	CUDASafeCallWithCleanUp(cudaMemcpy(g_pfPFBCoeff_d, g_pfPFBCoeff, sizePFB, cudaMemcpyHostToDevice));

	// allocate memory for FFT in and out arrays
	(void) fprintf(stdout, "\tAllocate memory for FFT arrays...\n");
	int sizeDataBlock = g_iNumSubBands * g_iNFFT * sizeof(float2);
	CUDASafeCallWithCleanUp(cudaMalloc((void **) &g_pf2FFTIn_d, sizeDataBlock));
	CUDASafeCallWithCleanUp(cudaMalloc((void **) &g_pf2FFTOut_d, sizeDataBlock));

	// set kernel parameters
	(void) fprintf(stdout, "\tSetting kernel parameters...\n");
	if(g_iNFFT < g_iMaxThreadsPerBlock) {
		g_dimBPFB.x  = g_iNFFT;
		g_dimBCopy.x = g_iNFFT;
	} else {
		g_dimBPFB.x  = g_iMaxThreadsPerBlock;
		g_dimBCopy.x = g_iMaxThreadsPerBlock;
	}
	g_dimGPFB.x  = (g_iNumSubBands * g_iNFFT) / g_dimBPFB.x;
	g_dimGCopy.x = (g_iNumSubBands * g_iNFFT) / g_dimBCopy.x;

	(void) fprintf(stdout, "\t\tKernel Parmaters are:\n\t\tgridDim(%d,%d,%d) blockDim(%d,%d,%d)\n",
							g_dimGPFB.x, g_dimGPFB.y, g_dimGPFB.z,
							g_dimBPFB.x, g_dimBPFB.y, g_dimGPFB.z);

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

	fprintf(stdout, "\nDevice for PFB successfully initialized!\n");
	return EXIT_SUCCESS;

}

int resetDevice() {
	cudaError_t cuErr = cudaDeviceReset();
	if (cuErr != cudaSuccess) {
		fprintf(stderr, "Device Reset Failed.\n");

		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

__global__ void map(char* dataIn,
			   		char2* dataOut,
			   		int channelSelect) 
{

	// select the channel range
	int channelMin = PFB_CHANNELS*channelSelect;
	
	int absIdx = 2 * blockDim.y*(blockIdx.x*DEF_NUM_CHANNELS + (channelMin+blockIdx.y)) + 2 * threadIdx.y;  // times 2 because we are mapping a sequence of values to char2 array.
	int mapIdx = blockDim.y*(blockIdx.x*gridDim.y + blockIdx.y) + threadIdx.y;

	dataOut[mapIdx].x = dataIn[absIdx];
	dataOut[mapIdx].y = dataIn[absIdx+1];
	return;
}

/* prepare data for PFB */
__global__ void PFB_kernel(char2* pc2Data,
                      float2* pf2FFTIn,
                      float* pfPFBCoeff)
{
    int i = (blockIdx.x * blockDim.x) + threadIdx.x;
    int iNFFT = (gridDim.x * blockDim.x);
    int j = 0;
    int iAbsIdx = 0;
    float2 f2PFBOut = make_float2(0.0, 0.0);
    char2 c2Data = make_char2(0, 0);

    for (j = 0; j < NUM_TAPS; ++j)
    {
        /* calculate the absolute index */
        iAbsIdx = (j * iNFFT) + i;

        /* get the address of the block */
        c2Data = pc2Data[iAbsIdx];
        
        f2PFBOut.x += (float) c2Data.x * pfPFBCoeff[iAbsIdx];
        f2PFBOut.y += (float) c2Data.y * pfPFBCoeff[iAbsIdx];
    }

    pf2FFTIn[i] = f2PFBOut;

    return;
}

__global__ void CopyDataForFFT(char2 *pc2Data, float2 *pf2FFTIn)
{
    int i = (blockIdx.x * blockDim.x) + threadIdx.x;

    pf2FFTIn[i].x = (float) pc2Data[i].x;
    pf2FFTIn[i].y = (float) pc2Data[i].y;

    return;
}

/* do fft on pfb data */
int doFFT()
{
    cufftResult iCUFFTRet = CUFFT_SUCCESS;

    /* execute plan */
    iCUFFTRet = cufftExecC2C(g_stPlan,
                             (cufftComplex*) g_pf2FFTIn_d,
                             (cufftComplex*) g_pf2FFTOut_d,
                             CUFFT_FORWARD);
    if (iCUFFTRet != CUFFT_SUCCESS)
    {
        (void) fprintf(stderr, "ERROR! FFT failed!\n");
        return EXIT_FAILURE;
    }

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
    if (g_pc2InBuf != NULL) {
        free(g_pc2InBuf);
        g_pc2InBuf = NULL;
    }
    if (g_pc2Data_d != NULL) {
        (void) cudaFree(g_pc2Data_d);
        g_pc2Data_d = NULL;
    }
    if (g_pf2FFTIn_d != NULL) {
        (void) cudaFree(g_pf2FFTIn_d);
        g_pf2FFTIn_d = NULL;
    }
    if (g_pf2FFTOut_d != NULL) {
        (void) cudaFree(g_pf2FFTOut_d);
        g_pf2FFTOut_d = NULL;
    }

    free(g_pfPFBCoeff);
    (void) cudaFree(g_pfPFBCoeff_d);

    /* destroy plan */
    /* TODO: check for plan */
    (void) cufftDestroy(g_stPlan);

    return;
}















