#include <stdio.h>
#include <math.h>
#include "../3rdparty/dSFMT-src-2.0/dSFMT.c"

double x1,x1p,x1preN;
double x2,x2p,x2preN;
double x3,x3p,x3preN;
double x4,x4p,x4preN;
double k1;
double k2;
double k3;
int x1ctr[100][5][5][5][5][5][5][5];
int x2ctr[100][5][5][5][5][5][5][5];
int x3ctr[100][5][5][5][5][5][5];
int x4ctr[100][5][5][5][5];

double fx1(double x1){
return -(k1*x1*x3-k2*x2)+k3*x2;
}
double fx2(double x2){
return (k1*x1*x3-k2*x2)-k3*x2;
}
double fx3(double x3){
return -(k1*x1*x3-k2*x2);
}
double fx4(double x4){
return k3*x2;
}


int discretize(double v, double xi[], int length){
for(int j=1;j<length-1;j++) if(v<xi[j]) return j-1;
return length-2;
}

int main(int argc, char *argv[]){
int myid=atoi(argv[1]);
double dt=0.01,halfdt=dt/2.0;
int tps=(int)(10.0/dt),t,i;
int block=tps/100,tb;
double halfF1,halfF2,F3,F4;
double x1i[]={0.0,3.0,6.0,9.0,12.0,15.0};
double x2i[]={0.0,3.0,6.0,9.0,12.0,15.0};
double x3i[]={0.0,3.0,6.0,9.0,12.0,15.0};
double x4i[]={0.0,3.0,6.0,9.0,12.0,15.0};

int x1pre,x1post,x1init=2;
int x2pre,x2post,x2init=0;
int x3pre,x3post,x3init=4;
int x4pre,x4post,x4init=0;
double k1i[]={0.0,0.2,0.4,0.6000000000000001,0.8,1.0};
double k2i[]={0.0,0.2,0.4,0.6000000000000001,0.8,1.0};
double k3i[]={0.0,0.2,0.4,0.6000000000000001,0.8,1.0};
int k1bin;
int k2bin;
int k3bin;



int sampleNo=1000;
dsfmt_t dsfmt;
int seed=7238+myid;
dsfmt_init_gen_rand(&dsfmt,seed);


for(int i=0;i<sampleNo;i++){
k1=0.0+dsfmt_genrand_close_open(&dsfmt)*1.0;
k2=0.0+dsfmt_genrand_close_open(&dsfmt)*1.0;
k3=0.0+dsfmt_genrand_close_open(&dsfmt)*1.0;
x1=x1i[x1init]+dsfmt_genrand_close_open(&dsfmt)*(x1i[x1init+1]-x1i[x1init]);
x2=x2i[x2init]+dsfmt_genrand_close_open(&dsfmt)*(x2i[x2init+1]-x2i[x2init]);
x3=x3i[x3init]+dsfmt_genrand_close_open(&dsfmt)*(x3i[x3init+1]-x3i[x3init]);
x4=x4i[x4init]+dsfmt_genrand_close_open(&dsfmt)*(x4i[x4init+1]-x4i[x4init]);
x1preN=x1;
x2preN=x2;
x3preN=x3;
x4preN=x4;
k1bin=discretize(k1,k1i,6);
k2bin=discretize(k2,k2i,6);
k3bin=discretize(k3,k3i,6);

for(int t=1;t<=tps;t++){
// x1
halfF1=halfdt*fx1(x1);
halfF2=halfdt*fx1(x1+halfF1);
F3=dt*fx1(x1+halfF2);
F4=dt*fx1(x1+F3);
x1p=x1+(2*halfF1+4*halfF2+2*F3+F4)/6.0;
// x2
halfF1=halfdt*fx2(x2);
halfF2=halfdt*fx2(x2+halfF1);
F3=dt*fx2(x2+halfF2);
F4=dt*fx2(x2+F3);
x2p=x2+(2*halfF1+4*halfF2+2*F3+F4)/6.0;
// x3
halfF1=halfdt*fx3(x3);
halfF2=halfdt*fx3(x3+halfF1);
F3=dt*fx3(x3+halfF2);
F4=dt*fx3(x3+F3);
x3p=x3+(2*halfF1+4*halfF2+2*F3+F4)/6.0;
// x4
halfF1=halfdt*fx4(x4);
halfF2=halfdt*fx4(x4+halfF1);
F3=dt*fx4(x4+halfF2);
F4=dt*fx4(x4+F3);
x4p=x4+(2*halfF1+4*halfF2+2*F3+F4)/6.0;
if(t%block==0){
tb=t/block-1;
x1pre=discretize(x1preN,x1i,6);
x2pre=discretize(x2preN,x2i,6);
x3pre=discretize(x3preN,x3i,6);
x4pre=discretize(x4preN,x4i,6);

x1post=discretize(x1,x1i,6);
x2post=discretize(x2,x2i,6);
x3post=discretize(x3,x3i,6);
x4post=discretize(x4,x4i,6);


x1ctr[tb][k1bin][k2bin][k3bin][x1pre][x2pre][x3pre][x1post]++;
x2ctr[tb][k1bin][k2bin][k3bin][x1pre][x2pre][x3pre][x2post]++;
x3ctr[tb][k1bin][k2bin][x1pre][x2pre][x3pre][x3post]++;
x4ctr[tb][k3bin][x2pre][x4pre][x4post]++;
x1preN=x1;
x2preN=x2;
x3preN=x3;
x4preN=x4;
}
x1=x1p;
x2=x2p;
x3=x3p;
x4=x4p;
}
}


// output
FILE *out;
char buffer[256];
snprintf(buffer, sizeof(buffer), "dummy.txt");
int idx=0;

for(tb=0;tb<100;tb++){
snprintf(buffer, sizeof(buffer), "../models/toy/batct/toyCTx1T%d_%d.txt", tb,myid);
out=fopen(buffer, "w");
idx=0;


for(int ki0=0;ki0<5;ki0++)
for(int ki1=0;ki1<5;ki1++)
for(int ki2=0;ki2<5;ki2++)
for(int vi0=0;vi0<5;vi0++)
for(int vi1=0;vi1<5;vi1++)
for(int vi2=0;vi2<5;vi2++)
for(int vi=0;vi<5;vi++)
{
int ctrtmp=(x1ctr[tb][ki0][ki1][ki2][vi0][vi1][vi2][vi]);
if(ctrtmp>0){
fprintf(out,"%d %d\n",idx,ctrtmp);
}
idx++;
}
fclose(out);
snprintf(buffer, sizeof(buffer), "../models/toy/batct/toyCTx2T%d_%d.txt", tb,myid);
out=fopen(buffer, "w");
idx=0;


for(int ki0=0;ki0<5;ki0++)
for(int ki1=0;ki1<5;ki1++)
for(int ki2=0;ki2<5;ki2++)
for(int vi0=0;vi0<5;vi0++)
for(int vi1=0;vi1<5;vi1++)
for(int vi2=0;vi2<5;vi2++)
for(int vi=0;vi<5;vi++)
{
int ctrtmp=(x2ctr[tb][ki0][ki1][ki2][vi0][vi1][vi2][vi]);
if(ctrtmp>0){
fprintf(out,"%d %d\n",idx,ctrtmp);
}
idx++;
}
fclose(out);
snprintf(buffer, sizeof(buffer), "../models/toy/batct/toyCTx3T%d_%d.txt", tb,myid);
out=fopen(buffer, "w");
idx=0;


for(int ki0=0;ki0<5;ki0++)
for(int ki1=0;ki1<5;ki1++)
for(int vi0=0;vi0<5;vi0++)
for(int vi1=0;vi1<5;vi1++)
for(int vi2=0;vi2<5;vi2++)
for(int vi=0;vi<5;vi++)
{
int ctrtmp=(x3ctr[tb][ki0][ki1][vi0][vi1][vi2][vi]);
if(ctrtmp>0){
fprintf(out,"%d %d\n",idx,ctrtmp);
}
idx++;
}
fclose(out);
snprintf(buffer, sizeof(buffer), "../models/toy/batct/toyCTx4T%d_%d.txt", tb,myid);
out=fopen(buffer, "w");
idx=0;


for(int ki0=0;ki0<5;ki0++)
for(int vi0=0;vi0<5;vi0++)
for(int vi1=0;vi1<5;vi1++)
for(int vi=0;vi<5;vi++)
{
int ctrtmp=(x4ctr[tb][ki0][vi0][vi1][vi]);
if(ctrtmp>0){
fprintf(out,"%d %d\n",idx,ctrtmp);
}
idx++;
}
fclose(out);
}

return 0;
}
