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

#include <float.h>		/* for FLT_MAX */

#define FALSE 					0
#define TRUE  					1
#define DEBUG					1
typedef unsigned char BYTE


#define DEF_CUDA_DEVICE			0

#define FILE_COEFF_PREFIX		"coeff"
#define FILE_COEFF_DATATYPE		"float"
#define FILE_COEFF_SUFFIX		".dat"

int loadDevice();
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