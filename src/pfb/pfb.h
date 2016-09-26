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

#define FALSE 					0
#define TRUE  					1
#define DEBUG					1

#define DEF_CUDA_DEVICE			0

#define DEF_LEN_SPEC			32
#define NUM_TAPS				8

// FFT Plan configuration
#define FFTPLAN_RANK 			1
#define FFTPLAN_ISTRIDE			(2*g_iNumSubBands)
#define FFTPLAN_OSTRIDE			(2*g_iNumSubBands)
#define FFTPLAN_IDIST			1
#define FFTPLAN_ODIST			1
#define FFPLAN_BATCH			(2 * g_iNumSubBands)

// coeff file configuration
#define FILE_COEFF_PREFIX		"coeff"
#define FILE_COEFF_DATATYPE		"float"
#define FILE_COEFF_SUFFIX		".dat"

typedef unsigned char BYTE;

// methods
int loadCoeff(int iCudaDevice);
int pfb();
int loadDataToMem(void);
int ReadData(void);

__global__ void doPFB();
__global__ void copyDataForFFT();

int doFFT();

void cleanUp(void);

#define CudaSafeCallWithCleanup(iRet) __CUDASafeCallWithCleanUp(iRet, __FILE__, __LINE__, &cleanUp)

void __CUDASafeCallWithCleanUp(cudaError_t iRet, const char* pcFile, const int iLine, void (*pcleanUp)(void));

#endif