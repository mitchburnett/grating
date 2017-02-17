/*
*   Kernels.h contains all cuda specific code and syntax.
*/

#include <cuda.h>
#include <cuda_runtime.h>

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