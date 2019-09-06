/*********************************************************************
 ** Stochastic Ranking Evolution Strategy                           **
 ** example 4 with SRES                                             **
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
 ** Date:   Mar 14, 2005;                                           **
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
  dim = 5;
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
  lb[0] = 78; ub[0] = 102;
  lb[1] = 33; ub[1] = 45;
  lb[2] = 27; ub[2] = 45;
  lb[3] = 27; ub[3] = 45;
  lb[4] = 27; ub[4] = 45;

  for(i=0;i<dim;i++)
    trsfm[i] = transform;

/********************************************************************
 ** end of user parameter setting                                  **
 ** get started                                                    **
 ********************************************************************/
  ESInitial(seed, &param, trsfm,fitness, es,constraint, dim, ub, lb,  \
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
  (*f) = 5.3578547*transform(x[2])*transform(x[2])+0.8356891*transform(x[0])*transform(x[4])+37.293239*transform(x[0])-40792.141;
  g[0] = 85.334407+0.0056858*transform(x[1])*transform(x[4])+0.0006262*transform(x[0])*transform(x[3])  \
         -0.0022053*transform(x[2])*transform(x[4])-92;;
  g[1] = -85.334407-0.0056858*transform(x[1])*transform(x[4])-0.0006262*transform(x[0])*transform(x[3])  \
         +0.0022053*transform(x[2])*transform(x[4]);
  g[2] = 80.51249+0.0071317*transform(x[1])*transform(x[4])+0.0029955*transform(x[0])*transform(x[1])  \
         +0.0021813*transform(x[2])*transform(x[2])-110;
  g[3] = -80.51249-0.0071317*transform(x[1])*transform(x[4])-0.0029955*transform(x[0])*transform(x[1])  \
         -0.0021813*transform(x[2])*transform(x[2])+90;
  g[4] = 9.300961+0.0047026*transform(x[2])*transform(x[4])+0.0012547*transform(x[0])*transform(x[2])  \
         +0.0019085*transform(x[2])*transform(x[3])-25;
  g[5] = -9.300961-0.0047026*transform(x[2])*transform(x[4])-0.0012547*transform(x[0])*transform(x[2])  \
         -0.0019085*transform(x[2])*transform(x[3])+20;

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

