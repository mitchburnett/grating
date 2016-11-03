#include "pfb.h"

__global__ map(char2 *dataIn,
			   char2 *dataOut,
			   int channelSelect) {

	/* map kernel
	The commented out section was when I was launching a thread for each index. The
	other code is where I launch only the threads I need and get the absIdx from that.
	i.e. The first code would require a kernel with the number of channels passed in as
	CHANNELS and the other requires only passing in PFB_CHANNELS.
	*/

	// select the channel range
	int channelMin = PFB_CHANNELS*channelSelect;
	//int channelMax = channelMin + (PFB_CHANNELS-1); //only need channel max when launching threads for CHANNELS otherwise the map works from channel min on up and gets the correct number of channels.

	/*
	// do noting if outside channels of interest
	int f = blockIdx.y;
	if ( f < channelMin || f > channelMax) {
		return;
	}
	// determine absolute index in dataIn
	f = f % PFB_CHANNELS;
	//int threadsPerBlock = blockDim.x*blockDim.y;
	int absIdx = blockDim.y*(blockIdx.x*gridDim.y + blockIdx.y) + threadIdx.y;
	int mapIdx = blockDim.y*(blockIdx.x*gridDim.y/PFB_CHANNELS + f) + threadIdx.y;
	*/
	int absIdx = blockDim.y*(blockIdx.x*CHANNELS + (channelMin+blockIdx.y)) + threadIdx.y;
	int mapIdx = blockDim.y*(blockIdx.x*gridDim.y + blockIdx.y) + threadIdx.y;

	dataOut[mapIdx] = dataIn[absIdx];
	return;

}
