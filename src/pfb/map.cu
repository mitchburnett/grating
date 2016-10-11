#include "pfb.h"

__global__ map(char2 *dataIn,
			   char2 *dataOut,
			   int channelSelect) {

	// select the channel range
	int channelMin = PFB_CHANNELS*channelSelect;
	int channelMax = channelMin + (PFB_CHANNELS-1);

	// do noting if outside channels of interest
	int f = blockIdx.y;
	if ( f < channelMin || f > channelMax) {
		return;
	}
	// determine absolute index in dataIn
	int threadsPerBlock = blockDim.x*blockDim.y;
	int absIdx = threadsPerBlock*(blockIdx.x*gridDim.y + blockIdx.y*blockDim.x) + threadIdx.y;
	int mapIdx = threadsPerBlock*(blockIdx.x*gridDim.y/PFB_CHANNELS + blockIdx.y*blockDim.x) + threadIdx.y;

	dataOut[mapIdx] = dataIn[absIdx];
	return;

}
