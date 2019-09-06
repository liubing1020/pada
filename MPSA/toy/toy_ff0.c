#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void load0(void);
void simulate(int kvar[]);
void output(void);
void outputMarginal(void);
double eval(int kvar[]);

int tps=100,tb,MAX_LINE=1024;
double dt=1.0;
double M[][5]={{-1,-1,-1,-1,-1,},
{1.0,3.0,5.0,7.0,9.0,},
{1.0,3.0,5.0,7.0,9.0,},
{1.5,4.5,7.5,10.5,13.5,},
{1.5,4.5,7.5,10.5,13.5,},
};

double x1prob[100][5][5][5][5][5][5][5];
double x2prob[100][5][5][5][5][5][5][5];
double x3prob[100][5][5][5][5][5][5];
double x4prob[100][5][5][5][5];
double result[5][101];
double x[5][101][5];

int P[]={1,2,3,};
int T[]={10,20,30,40,50,60,70,80,90,100,};
int X[]={1,2,3,4,};
double weights[]={0.642113,1.0,0.776737,0.435402,};
int k[]={-1,0,0,1,};
int xbinNum[]={-1,5,5,5,5,};
double expdata[10][4];

void load0(){

//loadExp
char buffer[MAX_LINE];
char *cp;
FILE *fin;

double tmp;
int i,j,k,ctr;
fin = fopen("../../SPA/models/toy/toydata.txt","r");
fgets(buffer, sizeof(buffer), fin); //header
fgets(buffer, sizeof(buffer), fin); //init
for(j=0;j<10;j++){
ctr=0;
fgets(buffer, sizeof(buffer), fin);
cp = (char *)strtok(buffer, "\t, "); //time
for(k=1;k<5;k++){
if(ctr==4) break;
cp = (char *)strtok(NULL, "\t, ");
tmp=strtod(cp, NULL);
if(k==X[ctr]){
expdata[j][ctr]=tmp;
ctr++;
}
}
}
fclose(fin);

// load CPD
int idx;double p;
int vb,vi0,vi1,vi2,vi3,vi4,vi5,vi6,vi7,vi8,vi9,vi10;
int pi0,pi1,pi2,pi3,pi4,pi5,pi6,pi7,pi8,pi9,pi10;
for(tb=0;tb<100;tb++){

snprintf(buffer, sizeof(buffer), "../../SPA/models/toy/ctn/tables/toyPx1T%d.txt", tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi2=idx%5;
idx=idx/5;
vi1=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
pi2=idx%5;
idx=idx/5;
pi1=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x1prob[tb][pi0][pi1][pi2][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "../../SPA/models/toy/ctn/tables/toyPx2T%d.txt", tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi2=idx%5;
idx=idx/5;
vi1=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
pi2=idx%5;
idx=idx/5;
pi1=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x2prob[tb][pi0][pi1][pi2][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "../../SPA/models/toy/ctn/tables/toyPx3T%d.txt", tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi2=idx%5;
idx=idx/5;
vi1=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
pi1=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x3prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "../../SPA/models/toy/ctn/tables/toyPx4T%d.txt", tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi1=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x4prob[tb][pi0][vi0][vi1][vb]=p;
}

fclose(fin);
}
}
void simulate(int kvar[]){
int t,b,tb,i,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15;
double sum=0,margin=1,norm=0;

for(i=0;i<3;i++) k[P[i]]=kvar[i];

//init
for(i=0;i<5;i++)
for(t=0;t<tps+1;t++)
for(b=0;b<5;b++)
x[i][t][b]=0;

for(i=0;i<5;i++)
for(t=0;t<tps+1;t++)
result[i][t]=0;


x[1][0][4]=1;
x[2][0][0]=1;
x[3][0][4]=1;
x[4][0][0]=1;

for(t=0;t<100;t++){
tb=t;
// normalization
for(i=1;i<5;i++){
norm=0; 
for(b=0;b<xbinNum[i];b++) norm+=x[i][t][b]; if(norm==0) norm=1;
for(b=0;b<xbinNum[i];b++) x[i][t][b]/=norm;
}


// normalization2
for(i=5;i<5;i++){
norm=0; 
for(b=0;b<xbinNum[i];b++) norm+=x[i][t][b];if(norm==0) norm=1;
for(b=0;b<xbinNum[i];b++) x[i][t][b]/=norm;
}

// x1
for(b=0;b<5;b++){
result[1][t]+=x[1][t][b]*M[1][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[1][t+1][b]+=x1prob[tb][k[1]][k[2]][k[3]][i0][i1][i2][b]*x[1][t][i0]*x[2][t][i1]*x[3][t][i2];
}
// x2
for(b=0;b<5;b++){
result[2][t]+=x[2][t][b]*M[2][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[2][t+1][b]+=x2prob[tb][k[1]][k[2]][k[3]][i0][i1][i2][b]*x[1][t][i0]*x[2][t][i1]*x[3][t][i2];
}
// x3
for(b=0;b<5;b++){
result[3][t]+=x[3][t][b]*M[3][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[3][t+1][b]+=x3prob[tb][k[1]][k[2]][i0][i1][i2][b]*x[1][t][i0]*x[2][t][i1]*x[3][t][i2];
}
// x4
for(b=0;b<5;b++){
result[4][t]+=x[4][t][b]*M[4][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
x[4][t+1][b]+=x4prob[tb][k[3]][i0][i1][b]*x[2][t][i0]*x[4][t][i1];
}
}
}

// output
void output(){
FILE *out;
char buffer[256];
snprintf(buffer, sizeof(buffer), "dummy.txt");
int idx=0,i,tb;

snprintf(buffer, sizeof(buffer), "toy_T.csv");
out=fopen(buffer, "w");
fprintf(out,"Time,");
for(i=1;i<5;i++)
fprintf(out,"x%d,",i);
fprintf(out,"\n");
for(tb=0;tb<tps;tb++){
fprintf(out,"%d,",tb);
for(i=1;i<5;i++)
fprintf(out,"%E,",result[i][tb]);
fprintf(out,"\n");
}
fclose(out);
}
// output marginal probabilities
void outputMarginal(){
FILE *out;
char buffer[256];
snprintf(buffer, sizeof(buffer), "dummy.txt");
int idx=0,t,i,b;

snprintf(buffer, sizeof(buffer), "toy_M.csv");
out=fopen(buffer, "w");
for(t=0;t<tps;t++)
for(i=1;i<5;i++)
for(b=0;b<xbinNum[i];b++)
fprintf(out,"%d,%d,%d,%E\n",t,i,b,x[i][t][b]);
fclose(out);
}
double eval(int kvar[]){
double err=0,value=0;int i,j;
simulate(kvar);
for(i=0;i<10;i++){
for(j=0;j<4;j++){
err=result[X[j]][T[i]]-expdata[i][j];
value+=weights[j]*err*err;
}
}
return value;
}
