#include "api.cuh"
#include "config.cu"
#include "stdio.h"

#define VAR(nr, min, max, bins, initbin) \
        if (threadIdx.y == (nr * WORK / NR_VARS)) {\
                localData[NR_CONSTS + 2 * NR_VARS + nr] = localData[NR_CONSTS + nr] = min  + ((real)initbin + extractNumber()) * (max - min) / (real)bins; \
                x_bins[nr + NR_VARS * threadIdx.x] = initbin;\
        }

#define AVAR(nr, min, max, bins)

#define CONST(nr, min, max, bins) \
        if (threadIdx.y == (nr * WORK / NR_CONSTS)) {\
                localData[nr] = extractNumber() * (max - min) + min;\
        }

#define EQUATION(rest...)
#define AEQUATION(rest...)

__device__ void GenRandom(unsigned int *x_bins, real *buffer_base, int thread_step) {
        real *localData = buffer_base;
        #include "model.cu"
}


#undef VAR
#undef AVAR
#undef CONST
#undef EQUATION
#undef AEQUATION

#define OPEN
#define CLOSE
#define C(c)
#define K(n)
#define X(n)
#define AX(n)
#define XX
#define CONST(rest...)
#define VAR(rest...)
#define AVAR(rest...)
#define ZERO

#define EQUATION(nr, formula, list) formula + 2
#define AEQUATION(rest...)

const int total_weight = 0
#include "model.cu"
;

#undef OPEN
#undef CLOSE
#undef C
#undef K
#undef X
#undef AX
#undef XX
#undef EQUATION
#undef AEQUATION
#undef CONST
#undef VAR
#undef AVAR
#undef ZERO

#define OPEN
#define CLOSE
#define C(c)
#define K(n)
#define X(n)
#define AX(n)
#define XX
#define CONST(rest...)
#define VAR(rest...)
#define AVAR(rest...)
#define ZERO

#define EQUATION(nr, formula, list) const int weight_ ## nr = 0 formula + 2;
#define AEQUATION(rest...)

#include "model.cu"

#undef MUL
#undef DIV
#undef PLUS
#undef MINUS
#undef OPEN
#undef CLOSE
#undef C
#undef K
#undef X
#undef AX
#undef XX
#undef EQUATION
#undef AEQUATION
#undef CONST
#undef VAR
#undef AVAR
#undef ZERO


#define HALFDT ( DT/ (real)2)

__device__ void FatNodes(real * buffer_base, int thread_step) {
	real __attribute__((unused)) * localData = buffer_base;

#define CONST(nr, rest...)
#define VAR(nr, rest...) 
#define AVAR(nr, rest...) 
#define C(c) (real)c
#define K(n) localData[n]
#define X(n) localData[NR_CONSTS + 2 * NR_VARS + n]
#define AX(n) localData[NR_CONSTS + 2 * NR_VARS + n]
#define XX xx_delta
#define DIV /
#define MUL *
#define PLUS +
#define MINUS -
#define OPEN (
#define CLOSE )
#define ZERO 0

#define AEQUATION(nr, formula, list) if (threadIdx.y  == (nr - NR_VARS) * WORK / NR_AVARS) localData[NR_CONSTS + 2 * NR_VARS + nr] = (formula); 

#define EQUATION(nr, formula, list)

	#include "model.cu"


#undef X
#undef AX
#undef C
#undef K
#undef CONST
#undef VAR
#undef AVAR
#undef EQUATION
#undef AEQUATION
#undef DIV
#undef MUL
#undef PLUS
#undef MINUS
#undef OPEN
#undef CLOSE
#undef ZERO

}

__shared__ volatile int sync_work[WORK * ((PAR_SIMS + 31) / 32)];
#define SYNC_WORKERS {if (threadIdx.x % 32 == 0) {sync_work[threadIdx.y * ((PAR_SIMS + 31) / 32) + threadIdx.x / 32]++; \
                        int syncs; \
                        do { \
                                for (syncs = (WORK - 1) * ((PAR_SIMS + 31) / 32); syncs >= 0; syncs--) \
                                        if (sync_work[syncs] + 1 == sync_work[threadIdx.y * ((PAR_SIMS + 31) / 32) + threadIdx.x / 32]) break; \
                        } while (syncs >= 0);}}



__device__ void RungeKutta(real * buffer_base, int thread_step, unsigned int * isOK) {
	real * localData = buffer_base;

#define CONST(nr, rest...)
#define VAR(nr, rest...)
#define AVAR(nr, rest...)
#define C(c) (real)c
#define K(n) localData[n]
#define X(n) localData[NR_CONSTS + n]
#define XX xx_delta
#define DIV /
#define MUL *
#define PLUS +
#define MINUS -
#define OPEN (
#define CLOSE )
#define ZERO 0

/*
The distribution to several threads is done by comparing the current_weight to the total_weight
	This weight is determined based on the operators weight, statically (use compiler macro redefinition to list the operators weight only)
*/

#define EQUATION(nr, formula, list) \
	if (threadIdx.y  == (current_weight * WORK / total_weight)) {\
		RK(localData[NR_CONSTS + NR_VARS + nr], X(nr), formula, localData[NR_CONSTS + 3*NR_VARS + nr], isOK[nr]); \
	}\
	current_weight = current_weight + weight_ ## nr;

#define AEQUATION(nr, formula, list)

	int current_weight = 0;
	#include "model.cu"


#undef X
#undef C
#undef K
#undef CONST
#undef VAR
#undef AVAR
#undef EQUATION
#undef AEQUATION
#undef DIV
#undef MUL
#undef PLUS
#undef MINUS
#undef OPEN
#undef CLOSE
#undef ZERO
}

#define CONST(rest...) 




#define VAR(nr, min, max, bins, initbin) if (threadIdx.y == (WORK + (nr * ACCESS / (NR_VARS + NR_AVARS)))) { \
	bin = 0; \
	real val = localData[NR_CONSTS + 2 * NR_VARS + nr];\
	for (int i = 1; i < bins; i++) { \
		if (val >= min + ((real)((max - min ) * i)) / (bins > 0 ? bins : 1)) bin = i;\
	} \
	if (bins > 0) newx_bins[nr + (NR_VARS + NR_AVARS) * threadIdx.x] = bin;\
}

#define AVAR(nr, min, max, bins) if (threadIdx.y == (WORK + (nr * ACCESS / (NR_VARS + NR_AVARS)))) { \
	bin = 0; \
	real val = localData[NR_CONSTS + 2 * NR_VARS + nr];\
	for (int i = 1; i < bins; i++) { \
		if (val >= min + ((real)((max - min) * i)) / (bins > 0 ? bins : 1)) bin = i;\
	} \
	if (bins > 0) newx_bins[nr + (NR_VARS + NR_AVARS) * threadIdx.x] = bin;\
}

#define EQUATION(rest...) 
#define AEQUATION(rest...)

__device__ void BinPackVars(unsigned int *newx_bins, real *buffer_base, int thread_step) {
	real * localData = buffer_base;
	int bin;
	#include "model.cu"
}
#undef CONST
#undef VAR
#undef AVAR

#define CONST(nr, min, max, bins) \
	bin = 0; \
	for (int i = 1; i < bins; i++) { \
		if (localData[nr] >= min + ((real)((max - min) * i)) / (bins > 0 ? bins : 1)) bin = i;\
	} \
	if (bins > 0) k_bins[nr + threadIdx.x * NR_CONSTS] = bin;
         

#define VAR(rest...) 
#define AVAR(rest...) 

__device__ void BinPackConsts(unsigned int *k_bins, real *buffer_base, int thread_step) {
	int bin;
	real * localData = buffer_base;
	#include "model.cu"
}

#undef CONST
#undef VAR
#undef AVAR
#undef EQUATION
#undef AEQUATION

#define CONST(rest...)
#define VAR(rest...)
#define AVAR(rest...)

#define K(n) [k_bins[n + NR_CONSTS * threadIdx.x]]

#define X(n) [x_bins[n + NR_VARS * threadIdx.x]]

#define AX(n) [newx_bins[n + (NR_VARS + NR_AVARS) * threadIdx.x]]

#undef offsetof
#define offsetof(st, m) ((size_t) ( (unsigned int *)&((st *)(0))->m - (unsigned int *)0 ))

#define EQUATION(nr, formula, list) if (valid && (threadIdx.y == WORK + (nr * ACCESS / (NR_VARS + NR_AVARS)))) incBin(data_out, offsetof(struct xctr, x ## nr ## ctr [block] list [newx_bins[nr + (NR_VARS + NR_AVARS) * threadIdx.x]]));

#define AEQUATION(nr, formula, list) if (valid && (threadIdx.y == WORK + (nr * ACCESS / (NR_VARS + NR_AVARS)))) incBin(data_out, offsetof(struct xctr, x ## nr ## ctr [block] list [newx_bins[nr + (NR_VARS + NR_AVARS) * threadIdx.x]]));


__device__ void IncBins(unsigned int *k_bins, unsigned int *x_bins, unsigned int *newx_bins, void * data_out, unsigned int block, int valid) {

	#include "model.cu"

}

#undef CONST
#undef VAR
#undef AVAR
#undef K
#undef X
#undef AX
#undef EQUATION
#undef AEQUATION

__device__ int VOTE( unsigned int * isOK)
{
	for(int i = 0; i< NR_VARS; i++)
		{
			if(isOK[i] == 0) return 0;
		}
	return 1;
}

__device__ void INIT_ST(real * buffer_base, real val)
{
	for(int i = 0; i< NR_VARS; i++)
	buffer_base[NR_CONSTS + 3 * NR_VARS + i] = (real)val;
}

__device__ void SET_VOTE(real * buffer_base, unsigned int * isOK)
{
	real voted_st = buffer_base[NR_CONSTS + 3 * NR_VARS];
	for(int i = 0; i< NR_VARS; i++)
	{
		if(voted_st > buffer_base[NR_CONSTS + 3 * NR_VARS + i])
			voted_st = buffer_base[NR_CONSTS + 3 * NR_VARS+ i];
	}
	buffer_base[NR_CONSTS + 4 * NR_VARS] = voted_st;
	INIT_ST(buffer_base, voted_st);
}



__global__ void __launch_bounds__(PAR_SIMS * (WORK + ACCESS), 1) kernel(void *data_out, int traces, int seed, int thread_step) {
	__shared__ unsigned int *k_bins, *x_bins, *newx_bins, *checkST;
	extern  __shared__ real shared_mem[];
	if (threadIdx.x == 0 && threadIdx.y == 0) {
		checkST = (unsigned int *)malloc(NR_VARS * PAR_SIMS * sizeof(unsigned int));
		k_bins = (unsigned int *)malloc(NR_CONSTS * PAR_SIMS * sizeof(unsigned int));
		x_bins = (unsigned int *)malloc(NR_VARS * PAR_SIMS * sizeof(unsigned int));
		newx_bins = (unsigned int *)malloc((NR_VARS + NR_AVARS) * PAR_SIMS * sizeof(unsigned int));
	}
	initializeGenerator((seed * gridDim.x + blockIdx.x) * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x);
	__syncthreads();
	
	for (int i = blockIdx.x * PAR_SIMS + threadIdx.x; (i / PAR_SIMS) * PAR_SIMS < traces; i += PAR_SIMS * gridDim.x) {
		real *buffer_base = shared_mem + threadIdx.x * thread_step;
		unsigned int *isOK = checkST + threadIdx.x * NR_VARS;
		if (threadIdx.y < WORK) {
			GenRandom(x_bins, buffer_base, thread_step);	
		}
			// assign X and K for all vars, also the bins for X
		__syncthreads();
		if (threadIdx.y == 0) {
			BinPackConsts(k_bins, buffer_base, thread_step);
		}
		__syncthreads();

	
		
		for (int block = -1; block < BLOCKS; block++) {
			if ((threadIdx.y < WORK) && (block < BLOCKS - 1) && (threadIdx.x % 32 == 0)) 
				sync_work[threadIdx.y * ((PAR_SIMS + 31) / 32) + threadIdx.x / 32] = 0;
			__syncthreads();
			if (threadIdx.y < WORK) {
				if (block < BLOCKS - 1) {
							real time = 0.0;
							real time_step = 1E-3;
							real INTERVAL = STEPS * ((real)(time_step));	
							INIT_ST(buffer_base,time_step);
							while(time < INTERVAL)	{
									int lim = 0;
									do{
										RungeKutta(buffer_base, thread_step,isOK);  					
										SYNC_WORKERS
										++ lim;
									} while ((!VOTE(isOK)) && (lim < 10));
									SET_VOTE(buffer_base,isOK);
									RungeKutta(buffer_base, thread_step,isOK); 
									time += buffer_base[NR_CONSTS + 4 * NR_VARS];
									SYNC_WORKERS
									INIT_ST(buffer_base,buffer_base[NR_CONSTS + 4 * NR_VARS]);
									for (int j = threadIdx.y; j < NR_VARS; j+= WORK) 
									buffer_base[NR_CONSTS + j] += buffer_base[NR_CONSTS + NR_VARS + j];
									SYNC_WORKERS
							}
					}
			} else {
				if (block >= 0) {
					BinPackVars(newx_bins, buffer_base, thread_step);
					IncBins(k_bins, x_bins, newx_bins, data_out, block, i < traces);
					for (int j = threadIdx.y - WORK; j < NR_VARS; j+= ACCESS) 
						x_bins[j + NR_VARS * threadIdx.x] = newx_bins[j + (NR_VARS + NR_AVARS) * threadIdx.x];
				}
			}
			if (block < BLOCKS - 1 && threadIdx.y < WORK) FatNodes(buffer_base, thread_step);

			__syncthreads();
			if (block < BLOCKS - 1) for (int j = threadIdx.y; j < NR_VARS; j+= WORK + ACCESS) 
				buffer_base[NR_CONSTS + 2 * NR_VARS + j] = buffer_base[NR_CONSTS + j];
			__syncthreads();
		}
	}
}



#define COALESCE(n) ((n) % 2 == 0 ? (n) + 1 : (n))


void simulation(struct cudaDeviceProp *props, void *data_out, int traces, unsigned int seed) {

	dim3 grid(props->multiProcessorCount);                                 // defines the number of blocks
        dim3 threads(PAR_SIMS, WORK + ACCESS);     // defines the pattern of threads inside each block
        int thread_step = COALESCE((NR_CONSTS + NR_VARS + NR_VARS + NR_VARS + NR_VARS + NR_AVARS + 1)); // coalescing is done at int level
        int sharedSize = PAR_SIMS * thread_step * sizeof(real);
        printf("Shared memory occupancy: %d (%.1f%%)\n", sharedSize, (float)sharedSize * 100 / props->sharedMemPerBlock);
		cudaFuncSetCacheConfig(kernel, cudaFuncCachePreferShared);
		cudaThreadSetLimit(cudaLimitMallocHeapSize, props->multiProcessorCount * 2 * PAR_SIMS * ((WORK + ACCESS) * 626 + NR_CONSTS + NR_VARS + NR_VARS + NR_VARS + NR_AVARS) * sizeof(unsigned int));
        kernel<<<grid, threads, sharedSize>>>(data_out, traces, seed, thread_step);
}


