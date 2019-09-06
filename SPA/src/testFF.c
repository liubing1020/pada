#include <stdio.h>
#include <stdlib.h>

#include "toy_ff0.c"

int main(int argc, char *argv[]){
  int i1,i2,i3;
  int k[3];
  
  // load CPTs
  load0();

  // for each combination of intervals of parameters, compute the error between 
  // data and prediction
  for(i1=0;i1<5;i1++)
    for(i2=0;i2<5;i2++)
      for(i3=0;i3<5;i3++){
	k[0]=i1;
	k[1]=i2;
	k[2]=i3;
	printf("%d%d%d %f\n",i1,i2,i3,eval(k));
      }
  

  // for a particular combination, output the mean and marginal after FF inference

  k[0]=0;
  k[1]=0;
  k[2]=1;
  simulate(k);

  output();
  outputMarginal();
  
  
}
