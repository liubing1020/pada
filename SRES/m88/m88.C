#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "sharefunc.h"
#include "ESSRSort.h"
#include "ESES.h"
#include "m88_ff_pe.c"

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
  dim = 163;
  miu = 2; //30
  lambda = 5; //200
  gen = 200;
  ub = NULL;
  lb = NULL;
  ub = ShareMallocM1d(dim);
  lb = ShareMallocM1d(dim);
  trsfm = (ESfcnTrsfm *)ShareMallocM1c(dim*sizeof(ESfcnTrsfm));
  //load0("/home/liubing/biomodel/bin/m88_centralized_sm_20_float_RK4/tmp/");
  load0("/home/liubing/biomodel/cluster/m88/tmp/");

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

  while(stats->curgen < param->gen)
    ESStep(population, param, stats, pf);


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
  int k[163];
  for(i=0;i<163;i++)
    k[i] = (int)x[i];
  
  value=obj(k);
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
