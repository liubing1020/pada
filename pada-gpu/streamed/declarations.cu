#include "api.cuh"
#include "config.cu"
#include "stdio.h"



#define HLFDT ( DT/ (real)2)

#define volatile 

__device__ struct struct_trace_data {
        real traceVars[NR_VARS];
        real traceConsts[NR_CONSTS];
        real traceFats[NR_AVARS + NR_VARS];
        real tracePost[NR_VARS];
        int k_bins[NR_CONSTS];
        int x_bins[NR_VARS];
        int newx_bins[NR_VARS + NR_AVARS];
};

__device__ void bin_pack(volatile struct struct_trace_data *this_trace, int nr, real val);


