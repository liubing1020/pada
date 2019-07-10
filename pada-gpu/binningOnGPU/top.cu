/**
   The top file of the simulation, loads the kernels after determining the runtime parameters
   The most important is derived as thread_step based on the number of constants and variables
 */

#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "multigpu.h"

// this is the prototype of a single GPU simulation
// the inner code does not make any assumptions on the data_out format, it carries the pointer onwards and passes it to the binning function
void simulation(struct cudaDeviceProp *props, void *data_out, int traces, unsigned int seed);

int main(int argc, char **argv) {
	int fileid;
	if (argc < 2) {
		printf("%s\n", "Error: At least the number of traces needs to be specified");
		return -1;
	}
	int TRACES = atoi(argv[1]);
	if (TRACES < 1) {
		printf("%s\n", "Error: Invalid number of traces");
		return -1;
	}
	int WHICH_DEVICE = 0;
	if (argc > 2) { 
		if (strcmp(argv[2], "-gpu") == 0) {
			if (argc < 4) {
				printf("%s\n", "Error: Invalid GPU specified");
				return -1;
			}
			WHICH_DEVICE = atoi(argv[3]);
			if (cudaSetDevice(WHICH_DEVICE) != cudaSuccess) {
				printf("%s\n", "Error: Invalid GPU specified");
				return -1;
			}
			if (argc>=5){ //./exec 3000 -gpu 0 -master 4
				useMultiGPU=1;
				if(strcmp(argv[4], "-master") == 0) {
					MSstatus=MASTER;
					useMultiGPU=atoi(argv[5]);
					if(useMultiGPU > MAX_GPU){
						printf("Error: at most %d GPUs available\n", MAX_GPU);
						return -1;
					}
				}
				else if (strcmp(argv[4], "-slave")==0){
					MSstatus=SLAVE;
					fileid=atoi(argv[5]);
				}
			}
		}
	}

	struct xctr *bins_host;
	struct xctr *bins_gpu;
	cudaEvent_t start, stop;
	bins_host = (struct xctr *)malloc(sizeof(struct xctr));
	if(!useMultiGPU || MSstatus==SLAVE){
		struct cudaDeviceProp props;	
		cudaGetDeviceProperties(&props, WHICH_DEVICE);

		cudaMalloc((void **)&bins_gpu, sizeof(struct xctr));
		cudaMemset(bins_gpu, 0, sizeof(struct xctr));

		cudaThreadSynchronize();
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start, 0);

		simulation(&props, bins_gpu, TRACES, WHICH_DEVICE);

		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		printf("Status: %s\n", cudaGetErrorString(cudaGetLastError()));
		cudaMemcpy(bins_host, bins_gpu, sizeof(struct xctr), cudaMemcpyDeviceToHost);
		cudaThreadSynchronize();
		float elapsedTime;
		cudaEventElapsedTime(&elapsedTime, start, stop);
		printf("Time: %.2fs\n", (double)elapsedTime/1000);
		cudaEventDestroy(start); 
		cudaEventDestroy(stop);
		cudaFree(bins_gpu);
	}

	//master thread will launch slave threads
	if(MSstatus==MASTER){
		printf("Using multi GPUs:\n");
		pthread_attr_t  attributes;
		pthread_attr_init(&attributes);
		pthread_attr_setdetachstate(&attributes, PTHREAD_CREATE_JOINABLE);

		pthread_t thread[MAX_GPU];
		ThreadData slavethd[MAX_GPU];
		ThreadData_master masterthd;

		//check gpus config
		if(gpudeviceinfo[0].device!=WHICH_DEVICE){ 
			printf("Please update multigpu.h\n");
			return -1;
		}
		masterthd._id=WHICH_DEVICE;
		masterthd._traces=TRACES/useMultiGPU;
		masterthd._data_h=bins_host;
		pthread_create(&thread[0], NULL, launchGPU_master, &masterthd);

		printf("master: device %d @ %s\n"
				,gpudeviceinfo[0].device, gpudeviceinfo[0].server);
		//create threads for slaves
		for(int i=1;i<useMultiGPU;i++){
			printf("slave: device %d @ %s\n"
					,gpudeviceinfo[i].device, gpudeviceinfo[i].server);
			slavethd[i]._id=i;
			slavethd[i]._MSstatus=SLAVE;
			//compare with master's server to check local/remote
			if(strcmp(gpudeviceinfo[i].server,gpudeviceinfo[0].server)==0){
				sprintf(slavethd[i]._cmd, "%s %d -gpu %d -slave %d"
						,argv[0], TRACES/useMultiGPU, gpudeviceinfo[i].device, i);
			}
			else{
				sprintf(slavethd[i]._cmd, "rsync -av %s %s@%s:%s/; ssh %s@%s 'cd %s; %s %d -gpu %d -slave %d'"
						, argv[0], USER, gpudeviceinfo[i].server, gpudeviceinfo[i].path
						, USER, gpudeviceinfo[i].server, gpudeviceinfo[i].path
						, argv[0], TRACES/useMultiGPU, gpudeviceinfo[i].device, i);
				sprintf(slavethd[i]._filecmd, "scp %s@%s:%s/bintemp%d ./"
						, USER, gpudeviceinfo[i].server, gpudeviceinfo[i].path, i);
			}
			pthread_create(&thread[i], NULL, launchGPU_slave, &slavethd[i]);
		}
		for(int i=0;i<useMultiGPU;i++){
			pthread_join(thread[i], NULL);
		}

		//collect output from slave gpus
		FILE *binpartfile;
		char binfilename[255];
		struct xctr *binpart;
		binpart=(struct xctr *)malloc(sizeof(struct xctr));
		//master gpu id is 0, the rest:     1,2,3 ...
		for(int i=1;i<useMultiGPU;i++){
			sprintf(binfilename,"bintemp%d",i);
			binpartfile=fopen(binfilename,"r");
			if(binpartfile!=NULL){
				fread(binpart,sizeof(struct xctr),1,binpartfile);
				fclose(binpartfile);
				unsigned int *u_part = ( unsigned int *)binpart;
				unsigned int *u_bin_host = ( unsigned int *)bins_host;
				for(int j=0;j<(sizeof(struct xctr)/sizeof(unsigned int));j++){
					u_bin_host[j]+=u_part[j];
				}
			}
		}
	}

	if(MSstatus==MASTER || !useMultiGPU){
		printf("Writing to files...\n");

		long checkup_counter = 0;

#include "WriteFiles.c"
		if (checkup_counter != (long)(NR_VARS + NR_AVARS) * BLOCKS * TRACES) printf("Checkup error on number of binned vars: %ld, should be %ld\n", checkup_counter, (long) (NR_VARS +NR_AVARS) * BLOCKS * TRACES);
		printf("Done\n");
	}
	else if (MSstatus==SLAVE){
		//slave gpu will write bin_host to a temp file
		FILE *tempfile;
		char tempfilename[255];
		sprintf(tempfilename, "bintemp%d", fileid);
		tempfile=fopen(tempfilename,"w");
		fwrite(bins_host,sizeof(struct xctr),1,tempfile);
		fclose(tempfile);
	}
	free(bins_host);

	return 0;
}
