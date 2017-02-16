#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <cuda.h>
#include <cuda_runtime.h> // Needed to include this or would complain about device identifier and cudaError_t and other cuda terms. (The question is this the correct file to include or are these terms included in cuda.h and I am just missing something to create the best wrapper.)

#define checkCudaErrors(err) __checkCudaErrors(err, #err, __FILE__, __LINE__)
void __checkCudaErrors(cudaError_t err, const char* const func, const char* file, const int line);

template<class T>
struct SharedMemory
{
    __device__ inline operator       T *()
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }

    __device__ inline operator const T *() const
    {
        extern __shared__ int __smem[];
        return (T *)__smem;
    }
};

__global__ void reduction(int* signal_d, int n);