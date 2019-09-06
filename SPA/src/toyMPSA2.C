#include <stdio.h>
#include <math.h>
#include "../../SPA/3rdparty/dSFMT-src-2.0/dSFMT.c"
/* header files and macros for CVODE */
#include "llnltyps.h"
#include "cvode.h"
#include "cvdense.h"
#include "nvector.h"
#include "dense.h"
#define Ith(v,i) N_VIth(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

int NEQ;
double RTOL, ATOL;
double T0, T1, Tm;
double x1;
double x2;
double x3;
double x4;
double k1;
double k2;
double k3;
double result[4][4];
double data[4][10];
int tn,xn;


double *P[]={&k1,&k2,&k3,};
int T[]={10,20,30,40,50,60,70,80,90,100,};
int X[]={1,2,3,4,};


static void difeq(integer N, real t, N_Vector y, N_Vector ydot, void *f_data);
int main(int argc, char *argv[]){
int ktest;
if(argc>1) ktest=atoi(argv[1]); else {printf("input ktest!\n"); exit(0);}
double lb[3];
double ub[3];
lb[0]=0.0; ub[0]=1.0;
lb[1]=0.0; ub[1]=1.0;
lb[2]=0.0; ub[2]=1.0;
/* set ODE parameters */
NEQ = 4;
RTOL = 1e-6;
ATOL = 1e-6;
T0 = 0.0;
Tm = 100;
tn=10;
xn=4;


x1=9.0;
x2=1.0;
x3=13.5;
x4=1.5;

k1=0.1;
k2=0.1;
k3=0.3;


// load data
char *cp;
FILE *in;
char buffer[256];
snprintf(buffer, sizeof(buffer), "../../SPA/models/toy/toydata.txt");
in = fopen(buffer, "r");
fgets(buffer, sizeof(buffer), in); //header
for(int j=0;j<10;j++){
int ctr2=0;
fgets(buffer, sizeof(buffer), in);
cp = (char *)strtok(buffer, "\t"); // time
for(int k=0;k<xn;k++){
cp = (char *)strtok(NULL, "\t");
double tmp = strtod(cp, NULL);
data[j][k]=tmp;
}
}


// output setting
FILE *out;
snprintf(buffer, sizeof(buffer), "Nk%d.txt",ktest);
out=fopen(buffer, "w");




int sampleNo=10000;
dsfmt_t dsfmt;
int seed=5589;
dsfmt_init_gen_rand(&dsfmt,seed);

double sumsum=0,threshold=200.0,freq=0;
int colNo=3;
int ctr[colNo];
for(int i=0;i<colNo;i++) ctr[i]=0;



for(int b=0;b<colNo;b++){
int validsample=0;
for(int i=0;i<sampleNo;i++){
k1=lb[0]+dsfmt_genrand_close_open(&dsfmt)*(ub[0]-lb[0]);
k2=lb[1]+dsfmt_genrand_close_open(&dsfmt)*(ub[1]-lb[1]);
k3=lb[2]+dsfmt_genrand_close_open(&dsfmt)*(ub[2]-lb[2]);
(*P[ktest])=lb[ktest]+(b+0.5)*(ub[ktest]-lb[ktest])/(colNo*1.0);


real ropt[OPT_SIZE], reltol, t, tout, ttmp;
long int iopt[OPT_SIZE];
N_Vector y;
real abstol;
void *cvode_mem;
int iout, flag, fail=0;
int i,k,j,l;
int t1;
double value = 0.0;
double sum,err,max,wmax;
double w[xn],m[xn];


y=N_VNew(NEQ, NULL);
Ith(y,1)=x1;
Ith(y,2)=x2;
Ith(y,3)=x3;
Ith(y,4)=x4;

reltol = RTOL;
abstol = ATOL;

cvode_mem = CVodeMalloc(NEQ, difeq, T0, y, BDF, NEWTON, SS, &reltol, &abstol, NULL, NULL, FALSE, iopt, ropt, NULL);

if(cvode_mem == NULL){
printf("CVodeMalloc failed.\n");
exit(1);
}
CVDense(cvode_mem, NULL, NULL);

ttmp=100;
for(iout=0;iout<tn;iout++){
tout=T[iout];
for(;ttmp<=tout;ttmp+=100){
flag = CVode(cvode_mem, ttmp, y, &t, NORMAL);
if (flag != SUCCESS){
printf("CVode failed, flag=%d.\n", flag);
fail=1;break;
}
}
if(fail==1) break;
result[0][iout] = Ith(y,1);
result[1][iout] = Ith(y,2);
result[2][iout] = Ith(y,3);
result[3][iout] = Ith(y,4);
ttmp=tout+100;
}
N_VFree(y);
CVodeFree(cvode_mem);
if(fail==0){
validsample++;
// evaluate
t1 = tn;
for(l=0,wmax=-1;l<xn;l++){
for(k=0,sum=0.0,max=-1; k<t1; k++){
sum+=data[l][k];
}
w[l]=t1/sum;
if(w[l]>wmax) wmax=w[l];
}

for(l=0;l<xn;l++){
for(k=0,sum=0.0; k<t1; k++){
err=data[l][k]-result[l][k];
sum +=err*err;
}
value+=sum*w[l]/wmax;
}

if(value<threshold){
ctr[b]++;
}
}
}
freq=ctr[b]/(validsample*1.0);
fprintf(out,"%f\t",freq);
}
fprintf(out,"\n");

fclose(out);
return 0;
}
/* differential equations */
static void difeq(integer N, real t, N_Vector y, N_Vector ydot, void *f_data){
real x1, dx1;
real x2, dx2;
real x3, dx3;
real x4, dx4;
x1 = Ith(y,1);
x2 = Ith(y,2);
x3 = Ith(y,3);
x4 = Ith(y,4);
dx1=-(k1*x1*x3-k2*x2)+k3*x2;
dx2=(k1*x1*x3-k2*x2)-k3*x2;
dx3=-(k1*x1*x3-k2*x2);
dx4=k3*x2;
Ith(ydot,1)=dx1;
Ith(ydot,2)=dx2;
Ith(ydot,3)=dx3;
Ith(ydot,4)=dx4;

return;
}

