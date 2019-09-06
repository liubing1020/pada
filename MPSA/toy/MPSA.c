/*
   input format: [par_id] [par_no] [dim] [threshold]
   e.g. Toy model: 
     $./MPSA 0 3 5 200.0 100
	 $./MPSA 1 3 5 200.0 100
	 $./MPSA 2 3 5 200.0 100
   output files: Nk0.txt, Nk1.txt, Nk2.txt
   these file are the input files for ksall.java
   note that the FF file must be included.
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "../../SPA/3rdparty/dSFMT-src-2.0/dSFMT.c"
#include "toy_ff0.c"

int main(int argc, char *argv[]){
  if(argc<6){
    printf("format: ./MPSA [par_id] [par_no] [dim] [threshold] [sampleNo]\n"); exit(0);
  }	
  int ktest=atoi(argv[1]);  // the parameter ID 
  int parNo=atoi(argv[2]);  // the number of unknow parameters
  int colNo=atoi(argv[3]); //  the number of intervals of unknown parameters 
  double threshold=strtod(argv[4], NULL); // the objective value threshold 200.0
  int sampleNo=atoi(argv[5]);
  int ub[parNo];
  for(int p=0;p<parNo;p++) ub[p]=colNo-1;
  int k[parNo];
  
  // load CPTs
  load0();
  
  // output setting
  FILE *out;
  char buffer[256];
  snprintf(buffer, sizeof(buffer), "Nk%d.txt",ktest);
  out=fopen(buffer, "w");

  dsfmt_t dsfmt;
  int seed=1234;
  dsfmt_init_gen_rand(&dsfmt,seed);
  
  double sumsum=0,freq=0;
  int ctr[colNo];
  for(int i=0;i<colNo;i++) ctr[i]=0;

  for(int i=1;i<=sampleNo;i++){
    for(int b=0;b<colNo;b++){
	  for(int p=0;p<parNo;p++)
	    k[p]=lround(dsfmt_genrand_close_open(&dsfmt)*ub[p]);      	  
	  k[ktest]=lround((b+0.5)*ub[ktest]/(colNo*1.0));

	  if(eval(k)<threshold){
	    ctr[b]++;
	  }
	  freq=ctr[b]/(i*1.0);
	  fprintf(out,"%f\t",freq);
	}
    fprintf(out,"\n");
  }
  fclose(out);
  return 0;
}
