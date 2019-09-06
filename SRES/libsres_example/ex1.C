/*********************************************************************
 ** Stochastic Ranking Evolution Strategy                           **
 ** example 1 with SRES                                             **
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
 ** Date:   Mar 7, 2005; Mar 8, 2005;                               **
 **                                                                 **
 ** Organization: Oak Ridge National Laboratory                     **
 ** Reference:                                                      **
 **   Thomas P. Runarsson and Xin Yao. 2000. Stochastic Ranking     **
 **   for Constrained Evolutionary Optimization. 4(3):284-294.      **
 **   http://cerium.raunvis.hi.is/~tpr/software/sres/               **
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>

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
  constraint = 9;
  dim = 13;
  miu = 30;
  lambda = 200;
  gen = 2000;
  ub = NULL;
  lb = NULL;
  ub = ShareMallocM1d(dim);
  lb = ShareMallocM1d(dim);
  trsfm = (ESfcnTrsfm *)ShareMallocM1c(dim*sizeof(ESfcnTrsfm));

/*********************************************************************
 ** set the up and low bounds on x here                             **
 ** lb[dim] and ub[dim]                                             **
 *********************************************************************/
  lb[0] = 0; ub[0] = 1;
  lb[1] = 0; ub[1] = 1;
  lb[2] = 0; ub[2] = 1;
  lb[3] = 0; ub[3] = 1;
  lb[4] = 0; ub[4] = 1;
  lb[5] = 0; ub[5] = 1;
  lb[6] = 0; ub[6] = 1;
  lb[7] = 0; ub[7] = 1;
  lb[8] = 0; ub[8] = 1;
  lb[9] = 0; ub[9] = 100;
  lb[10] = 0; ub[10] = 100;
  lb[11] = 0; ub[11] = 100;
  lb[12] = 0; ub[12] = 1;

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

  (*f) = 5*(transform(x[0])+transform(x[1])+transform(x[2])+transform(x[3]))   \
         - 5*(transform(x[0])*transform(x[0])+transform(x[1])*transform(x[1])  \
         +transform(x[2])*transform(x[2])+transform(x[3])*transform(x[3]))  \
         - (transform(x[4])+transform(x[5])+transform(x[6])+transform(x[7])  \
         +transform(x[8])+transform(x[9])+transform(x[10])+transform(x[11])  \
         +transform(x[12]));
  g[0] = 2*transform(x[0])+2*transform(x[1])+transform(x[9])+transform(x[10])-10;
  g[1] = 2*transform(x[0])+2*transform(x[2])+transform(x[9])+transform(x[11])-10;
  g[2] = 2*transform(x[1])+2*transform(x[2])+transform(x[10])+transform(x[11])-10;
  g[3] = -8*transform(x[0])+transform(x[9]);
  g[4] = -8*transform(x[1])+transform(x[10]);
  g[5] = -8*transform(x[2])+transform(x[11]);
  g[6] = -2*transform(x[3])-transform(x[4])+transform(x[9]);
  g[7] = -2*transform(x[5])-transform(x[6])+transform(x[10]);
  g[8] = -2*transform(x[7])-transform(x[8])+transform(x[11]);

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
