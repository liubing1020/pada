// receives the offset in the table where the increment should take place //
// thread_id = [0, BIN_THREADS] is guaranteed to be different for all the threads that call incBin in parallel

__device__ void incBin(void *data_out, unsigned int offset, int thread_id) {
	unsigned int *where = ((unsigned int *)data_out) + offset;
	atomicAdd(where, 1);
}
