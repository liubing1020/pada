#include "binning_structure.h"
#include "rng/GenRandom.cu"
#include "rk.cu"

__device__ void incBin(void *data_out, unsigned int offset, int thread_id);
