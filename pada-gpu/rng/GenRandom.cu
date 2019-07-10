/*
(c) Andrei Hagiescu 2010
This file generates 'real' random numbers in a specified interval
It uses a Mersene RNG whose source code is specified on Wikipedia
*/

__shared__ unsigned int *indx, *MT;

__device__ void initializeGenerator(unsigned int which) {
	int threadId = threadIdx.y * blockDim.x + threadIdx.x;
	if (threadId == 0) {
		indx = (unsigned int *)malloc(blockDim.x * blockDim.y * sizeof(unsigned int));
		MT = (unsigned int *)malloc(blockDim.x * blockDim.y * 625 * sizeof(unsigned int));
	}
	__syncthreads();
	MT[threadId * 625] = which;
	for (int i = 1; i < 623; i++) MT[threadId * 625 + i] = 0x6c078965 * (MT[threadId * 625 + i-1] ^ (MT[threadId * 625 + i-1] >> 30)) + i;
	indx[threadId] = 0;
}
 

 // Generate an array of 624 untempered numbers
__device__ void generateNumbers() {
	int threadId = threadIdx.y * blockDim.x + threadIdx.x;
	for (int i = 0; i < 624; i++) {
		unsigned int y = (MT[threadId * 625 + i] & 0x01) + (MT[threadId * 625 + (i+1) % 624] & 0x7FFFFFFF);
		MT[threadId * 625 + i] = MT[threadId * 625 + (i + 397) % 624] ^ (y >> 1);
		if (y % 2 == 1) { // y is odd
			MT[threadId * 625 + i] = MT[threadId * 625 + i] ^ 0x9908b0df;
		}
	}
}


__device__ real extractNumber() {
	int threadId = threadIdx.y * blockDim.x + threadIdx.x;
	if (indx[threadId] == 0) generateNumbers();

	unsigned int y = MT[threadId * 625 + indx[threadId]];
	y = y ^ (y >> 11);
	y = y ^ ((y << 7) & 0x9d2c5680);
	y = y ^ ((y << 15) & 0xefc60000);
	y = y ^ (y >> 18);
	indx[threadId]++;
	if (indx[threadId] == 624) indx[threadId] = 0;
	return (real)y / (real)0xffffffff;
}

