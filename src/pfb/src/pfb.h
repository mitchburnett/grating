#ifndef __PFB_H__
#define __PFB_H__

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cufft.h>

#include <string.h>		/* for strncopy(), memcpy(), strerror()*/
#include <sys/types.h>  /* for open()  */
#include <sys/stat.h>	/* for open()  */
#include <fcntl.h>		/* for open()  */
#include <unistd.h>		/* for close() */
#include <errno.h>		/* for errno   */

#include <float.h>		/* for FLT_MAX */

#include <getopt.h>		/* for option parsing */

#include <python2.7/Python.h> /* for executing coeff gen file */

#define FALSE 					0
#define TRUE  					1
#define DEBUG					1

#define DEF_CUDA_DEVICE			0

#define DEF_SIZE_READ			262144	// data block size. should this be set dynamically once I get the data?
#define DEF_LEN_SPEC			32	// Transform size
#define NUM_TAPS				8	// PFB Decimation factor
#define DEF_NUM_CHANNELS		25  // System spec for total number of channels
#define PFB_CHANNELS			5	// Number of coarse channels through PFB
#define DEF_NUM_ELEMENTS		64  // System spec for number of elements
#define SAMPLES					4000// Time samples.

// FFT Plan configuration
#define FFTPLAN_RANK 			1				 // dimension of the transform
#define FFTPLAN_ISTRIDE			(g_iNumSubBands) // The distance between two successive input time elements. - (polarization*numsubbands).
#define FFTPLAN_OSTRIDE			(g_iNumSubBands) // Similar to ostride to maintain data structure
#define FFTPLAN_IDIST			1				 // The distance between the first elements of two consecutive batches in the input data. Each FFT operation is a 'batch'. Each subband is a time series and we need a FFT for each subband. Since we have interleaved samples the distance between consecutive batches is 1 sample.
#define FFTPLAN_ODIST			1				 // Simailar to odist to maintian data structure
#define FFTPLAN_BATCH			(g_iNumSubBands) // The total number of FFTs to perform per call to DoFFT().

// coeff file configuration
#define FILE_COEFF_PREFIX		"coeff"
#define FILE_COEFF_DATATYPE		"float"
#define FILE_COEFF_SUFFIX		".dat"

typedef unsigned char BYTE;

struct params{
	int nfft;
	int taps;
	int subbands;
	int select;
};


// methods
int loadCoeff(int iCudaDevice);
int runPFB(char* inputData_h, float2* outputData_h, params pfbParams);
//int loadDataToMem(void);
//int ReadData(void);

__global__ void PFB_kernel(char2* pc2Data, float2* pf2FFTIn, float* pfPFBCoeff);
__global__ void map(char* dataIn, char2* dataOut, int channelSelect);
__global__ void CopyDataForFFT(char2* pc2Data, float2* pf2FFTIn);

int doFFT();

void cleanUp(void);
int resetDevice(void);

void genCoeff(int argc, char* argv[]);

#define CUDASafeCallWithCleanUp(iRet) __CUDASafeCallWithCleanUp(iRet, __FILE__, __LINE__, &cleanUp)

void __CUDASafeCallWithCleanUp(cudaError_t iRet, const char* pcFile, const int iLine, void (*pcleanUp)(void));

#endif