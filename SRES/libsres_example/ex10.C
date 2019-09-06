/*********************************************************************
 ** Stochastic Ranking Evolution Strategy                           **
 ** example 10 with SRES                                             **
 **                                                                 **
 ** For ACADEMIC RESEARCH, this is licensed with GPL license        **
 ** For COMMERCIAL ACTIVITIES, please contact the authors           **
 **                                                                 **
 ** Copyright (C) 2005 Xinglai Ji (jix1@ornl.gov)                   **
 **                                                                 **
 ** This program is distributed in the hope that it will be useful, **
 ** but WITHOUT ANY WARRANTY; without even the implied warranty of  **
 ** MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the    **
 ** GNU General Public License for more details.                    **
 **                                                                 **
 ** You should have received a copy of the GNU General Public       **
 ** License along with is program; if not, write to the Free        **
 ** Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, **
 ** MA 02111-1307, USA.                                             **
 **                                                                 **
 ** Author: Xinglai Ji (jix1@ornl.gov)                              **
 ** Date:   Mar 15, 2005;                                           **
 **                                                                 **
 ** Organization: Oak Ridge National Laboratory                     **
 ** Reference:                                                      **
 **   Thomas P. Runarsson and Xin Yao. 2000. Stochastic Ranking     **
 **   for Constrained Evolutionary Optimization. 4(3):284-294.      **
 **   http://cerium.raunvis.hi.is/~tpr/software/sres/               **
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sharefunc.h"
#include "ESSRSort.h"
#include "ESES.h"

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
  gamma = esDefGamma;
  alpha = esDefAlpha;
  varphi = esDefVarphi;
  retry = esDefRetry;

  pf = essrDefPf;

  es = esDefESSlash;
  constraint = 6;
  dim = 8;
  miu = 30;
  lambda = 200;
  gen = 1750;
  ub = NULL;
  lb = NULL;
  ub = ShareMallocM1d(dim);
  lb = ShareMallocM1d(dim);
  trsfm = (ESfcnTrsfm *)ShareMallocM1c(dim*sizeof(ESfcnTrsfm));

/*********************************************************************
 ** set the up and low bounds on x here                             **
 ** lb[dim] and ub[dim]                                             **
 *********************************************************************/
/*
  lb[0] = 100; ub[0] = 10000;
  lb[1] = 1000; ub[1] = 10000;
  lb[2] = 1000; ub[2] = 10000;
  for(i=3;i<dim;i++)
  {
    lb[i] = 10; 
    ub[i] = 1000;
  }
*/
  lb[0] = 2; ub[0] = 4;
  lb[1] = 3; ub[1] = 4;
  lb[2] = 3; ub[2] = 4;
  for(i=3;i<dim;i++)
  {
    lb[i] = 1; 
    ub[i] = 3;
  }

  for(i=0;i<dim;i++)
    trsfm[i] = transform;

/********************************************************************
 ** end of user parameter setting                                  **
 ** get started                                                    **
 ********************************************************************/
  ESInitial(seed, &param, trsfm,fitness,es, constraint, dim, ub, lb,  \
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
  (*f) = transform(x[0])+transform(x[1])+transform(x[2]);
  g[0] = -1+0.0025*(transform(x[3])+transform(x[5]));
  g[1] = -1+0.0025*(transform(x[4])+transform(x[6])-transform(x[3]));
  g[2] = -1+0.01*(transform(x[7])-transform(x[4]));
  g[3] = -transform(x[0])*transform(x[5])+833.33252*transform(x[3])+100*transform(x[0])-83333.333;
  g[4] = -transform(x[1])*transform(x[6])+1250*transform(x[4])+transform(x[1])*transform(x[3])-1250*transform(x[3]);
  g[5] = -transform(x[2])*transform(x[7])+1250000+transform(x[2])*transform(x[4])-2500*transform(x[4]);

  return;
}

/*********************************************************************
 ** double transform(double x)                                      **
 ** transform x: for coding                                         **
 *********************************************************************/
double transform(double x)
{
  double y;
                                                                                
  y = pow(10,x);
                                                                                
  return y;
}
