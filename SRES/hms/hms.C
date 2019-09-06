#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "sharefunc.h"
#include "ESSRSort.h"
#include "ESES.h"
#include "hms_ff_pe.c"

/*********************************************************************
 ** fitness(double *x, double *f, double *g)                        **
 ** change the function name as you want                            **
 ** x[dim], g[constraint]                                           **
 *********************************************************************/
void fitness(double *, double *, double *);

/*********************************************************************
 ** double transform(double x)                                      **
 ** transform x: for coding                                         **
 *********************************************************************/
double transform(double);

int main(int argc, char ** argv)
{
  FILE *out,*input2,*out2;
  char buffer[64];
  int i;
  ESParameter *param;
  ESPopulation *population;
  ESStatistics *stats;
  ESfcnTrsfm *trsfm;
  unsigned int seed;
  int es;
  int constraint, dim;
  double *ub, *lb;
  int miu, lambda, gen;
  double gamma, alpha, varphi;
  int retry;
  double pf;

/*********************************************************************
 ** change the parameters here as you want                          **
 ** random seed, gamma, alpha, varphi, retry number, pf,            **
 ** if you dont know how to set, keep default settings              **
 **                                                                 **
 ** constraint number, x dimension, miu:parent number,              **
 ** lambda:offspring number, generation number,                     **
 ** up and low bounds on x,                                         **
 *********************************************************************/
  seed = shareDefSeed;
  gamma = 0.85;//esDefGamma;
  alpha = 0.2;//esDefAlpha;
  varphi = 1; //esDefVarphi;
  retry = esDefRetry;

  pf = essrDefPf;

  es = esDefESSlash;
  constraint = 1;
  dim = 14;
  
  input2=fopen("input2.dat", "r");
  fgets(buffer, sizeof(buffer), input2);
  miu=strtol(buffer,NULL,0);
  fgets(buffer, sizeof(buffer), input2);
  lambda=strtol(buffer,NULL,0);
  fgets(buffer, sizeof(buffer), input2);
  gen=strtol(buffer,NULL,0);
  fclose(input2);  
  //miu = 30; //30
  //lambda = 200; //200
  //gen = 200;
  ub = NULL;
  lb = NULL;
  ub = ShareMallocM1d(dim);
  lb = ShareMallocM1d(dim);
  trsfm = (ESfcnTrsfm *)ShareMallocM1c(dim*sizeof(ESfcnTrsfm));
  load9("/hpctmp/dcsliubi/hms/tmp/");

/*********************************************************************
 ** set the up and low bounds on x here                             **
 ** lb[dim] and ub[dim]                                             **
 *********************************************************************/
  
  for(i=0;i<dim;i++){
    lb[i] = 0; ub[i] = 4;
  }



  for(i=0;i<dim;i++)
    trsfm[i] = transform;

/********************************************************************
 ** end of user parameter setting                                  **
 ** get started                                                    **
 ********************************************************************/
  ESInitial(seed, &param, trsfm, fitness, es,constraint, dim, ub, lb,  \
            miu, lambda, gen, gamma, alpha, varphi,  \
            retry, &population, &stats);

  out=fopen("output1.dat", "w");
  out2=fopen("output2.dat", "w");			

  while(stats->curgen < param->gen){
    ESStep(population, param, stats, pf);
    if(stats->curgen % 10 ==0)
      fprintf(out2,"gen=%d,dt=%d,bestgen=%d,bestfitness=%f,phi=%f,\n",stats->curgen,stats->dt,stats->bestgen,stats->bestindvdl->f,stats->bestindvdl->phi);
  }
  fprintf(out2,"gen=%d,dt=%d,bestgen=%d,bestfitness=%f,phi=%f,\n",stats->curgen,stats->dt,stats->bestgen,stats->bestindvdl->f,stats->bestindvdl->phi);
  for(i=0;i<dim;i++)
    fprintf(out,"%f,%f\n",stats->bestindvdl->op[i],stats->bestindvdl->sp[i]);
  fclose(out);
  fclose(out2);
  
  /*		
  while(stats->curgen < param->gen)
    ESStep(population, param, stats, pf);
*/

  ESDeInitial(param, population, stats);

  ShareFreeM1d(ub);
  ub = NULL;
  ShareFreeM1d(lb);
  lb = NULL;
  ShareFreeM1c((char*)trsfm);
  trsfm = NULL;

  return 0;
}


/*********************************************************************
 ** set the fitness function here                                   **
 ** x[dim], g[constraint]                                           **
 *********************************************************************/
void fitness(double *x, double *f, double *g)
{
  double value = 0.0;
  int i;
  int k[14];
  for(i=0;i<14;i++)
    k[i] = (int)x[i];

  value=eval(k);
  //value=obj(k);
  //printf("%d %d %d = %f\n",k[0],k[1],k[2],value);
  //printf("%E\n",value);

  (*f) = value;
  g[0] = 0.0;

  return;
}

/*********************************************************************
 ** double transform(double x)                                      **
 ** transform x: for coding                                         **
 *********************************************************************/
double transform(double x)
{
  double y;

  y = x;

  return y;
}
