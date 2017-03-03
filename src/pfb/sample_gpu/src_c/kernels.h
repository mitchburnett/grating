/*
*   Kernels.h contains all cuda specific code and syntax.
*/

#include <cuda.h>
#include <cuda_runtime.h>

// struct SharedMemory
// {
//     __device__ inline operator       int *()
//     {
//         extern __shared__ int __smem[];
//         return (int *)__smem;
//     }

//     __device__ inline operator const int *() const
//     {
//         extern __shared__ int __smem[];
//         return (int *)__smem;
//     }
// };

__global__ void reduction(int* signal_d, int n);