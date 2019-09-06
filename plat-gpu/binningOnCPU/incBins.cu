// receives the offset in the table where the increment should take place
// thread_id = [0, BIN_THREADS] is guaranteed to be different for all the threads that call incBin in parallel

__device__ void incBin(void *data_out, unsigned int offset, int thread_id) {
	unsigned int *CallCtr=(unsigned int*)data_out;
	unsigned int *u_data_out=(unsigned int*)data_out+BIN_THREADS;
	u_data_out[CallCtr[thread_id]*BIN_THREADS+thread_id]=offset+1;
	CallCtr[thread_id]++;
}
