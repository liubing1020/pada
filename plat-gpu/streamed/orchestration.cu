__device__ void bin_pack(volatile struct struct_trace_data *this_trace, int nr, real val) {
        int bin = 0;
        for (int i = 1; i < var_bins[nr]; i++) {
                if (val >= var_min[nr] + ((real)((var_max[nr] - var_min[nr]) * i)) / var_bins[nr]) bin = i;
        }
        this_trace -> newx_bins[nr] = bin;
}

__device__ real get_random(real min, real max) {
        return extractNumber()  * (max - min) + min;
}

/* main entry point */
/* launch bounds to assist the compiler regarding the maximum number of registers per thread */

__global__ __launch_bounds__(nrWarps * 32, 1) void kernel(void * data_out, int traces, unsigned int seed) {
	__shared__ struct_trace_data *trace_data;
        volatile extern __shared__ real sharedData[];
	if (threadIdx.x == 0 && threadIdx.y == 0) {
		trace_data = (struct_trace_data *)malloc(32 * sizeof(struct_trace_data));
	}
	initializeGenerator((seed * gridDim.x + blockIdx.x) * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x);
        __syncthreads();

        /*   Each trace requires to have its working set in global memory;   **
        **   it will be prefetched as needed during computation              */

        volatile struct struct_trace_data *this_trace = &trace_data[threadIdx.x];

        for (int trace = blockIdx.x * 32 + threadIdx.x; (trace / 32) * 32 < traces ; trace += 32 * gridDim.x) {
                int i;

                if (threadIdx.y == 0) for (i = 0; i < NR_CONSTS; i++) {
			int bin = 0;
                        real val = get_random(const_min[i], const_max[i]);
                        this_trace -> traceConsts[i] = val;
                        for (int j = 1; j < const_bins[i]; j++) {
                                if (val >= const_min[i] + ((real)((const_max[i] - const_min[i]) * j)) / const_bins[i]) bin = j;
                        }
                        this_trace -> k_bins[i] = bin;
                }

                if (threadIdx.y == 0) for (i = 0; i < NR_VARS; i++) {
                        this_trace -> x_bins[i] = var_initbin[i];
                        this_trace -> traceVars[i] =
                                get_random(var_min[i] + (var_max[i] - var_min[i]) * (real)var_initbin[i] / (real)var_bins[i],
                                           var_min[i] + (var_max[i] - var_min[i]) * (real)(var_initbin[i] + 1) / (real)var_bins[i]);
                }

                __syncthreads();
		__threadfence();

                int block, t;

                for (block = 0; block < BLOCKS; block++) {
                        for (t = 0; t < STEPS; t++) {
                                compute(t, this_trace, sharedData + threadIdx.x * THREAD_STEP);
                                __syncthreads();
				__threadfence();
				for (i = threadIdx.y; i < NR_VARS; i += nrWarps) this_trace -> traceVars[i] = this_trace -> tracePost[i];
                                __syncthreads();
				__threadfence();
                	}



#undef offsetof
#define offsetof(st, m) ((size_t) ( (unsigned int *)&((st *)(0))->m - (unsigned int *)0 ))

#define BIN_THREADS (32 * gridDim.x * nrWarps)
#include "incBins.cu"

                        #define VAR(...)
                        #define AVAR(...)
                        #define CONST(...)
                        #define K(n) [this_trace -> k_bins[n]]
                        #define X(n) [this_trace -> x_bins[n]]
                        #define AX(n) [this_trace -> newx_bins[n]]
                        #define EQUATION(nr, formula, list) \
                                if (trace < traces && (nr % nrWarps == threadIdx.y)) incBin(data_out, offsetof(struct xctr, x ## nr ## ctr [block] list [this_trace -> newx_bins[nr]]), (blockIdx.x * 32 + threadIdx.x) * nrWarps + threadIdx.y);


                        #define AEQUATION(nr, formula, list) \
                                if (trace < traces && (nr % nrWarps == threadIdx.y)) incBin(data_out, offsetof(struct xctr, x ## nr ## ctr [block] list [this_trace -> newx_bins[nr]]), (blockIdx.x * 32 + threadIdx.x) * nrWarps + threadIdx.y);

                        #include "model.cu"

                        #undef X
                        #undef AX
                        #undef K
                        #undef EQUATION
                        #undef AEQUATION
                        #undef VAR
                        #undef AVAR
                        #undef CONST

                        __syncthreads();
			__threadfence();

                        for (i = threadIdx.y; i < NR_VARS; i += nrWarps) this_trace -> x_bins[i] = this_trace -> newx_bins[i];

                        __syncthreads();
			__threadfence();
                }
        }
}


void simulation(struct cudaDeviceProp *props, void *data_out, int traces, unsigned int seed) {
	dim3 grid(props->multiProcessorCount);                                 // defines the number of blocks
        dim3 threads(32, nrWarps);     // defines the pattern of threads inside each block
	cudaFuncSetCacheConfig(kernel, cudaFuncCachePreferShared);
	if (cudaSuccess != cudaThreadSetLimit(cudaLimitMallocHeapSize, props->multiProcessorCount * 8 * 32 * (nrWarps * 626 * sizeof(unsigned int) + sizeof(struct_trace_data)))) printf("Error allocating memory");
	cudaThreadSynchronize();
        kernel<<<grid, threads, 32 * THREAD_STEP * sizeof(real)>>>(data_out, traces, seed);
}


