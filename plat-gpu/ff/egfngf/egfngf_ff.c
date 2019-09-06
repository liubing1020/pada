#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void load0(char*);
void simulate(int kval[]);
void output(void);
void outputMarginal(void);
double eval(int kval[]);

int tps=100,tb,MAX_LINE=1024;
float dt=1;
double M[][5]={{-1,-1,-1,-1,-1,},
{1100220.0,3300660.0,5501100.0,7701540.0,9901980.0,},
{50160.0,150480.0,250800.0,351120.0,451440.0,},
{8800.0,26400.0,44000.0,61600.0,79200.0,},
{8800.0,26400.0,44000.0,61600.0,79200.0,},
{1100.0,3300.0,5500.0,7700.0,9900.0,},
{1100.0,3300.0,5500.0,7700.0,9900.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{66000.0,198000.0,330000.0,462000.0,594000.0,},
{66000.0,198000.0,330000.0,462000.0,594000.0,},
{66000.0,198000.0,330000.0,462000.0,594000.0,},
{66000.0,198000.0,330000.0,462000.0,594000.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
{13200.0,39600.0,66000.0,92400.0,118800.0,},
};

double x1prob[100][5][5][5][5][5][5];
double x2prob[100][5][5][5][5][5][5];
double x3prob[100][5][5][5][5][5][5];
double x4prob[100][5][5][5][5][5][5];
double x5prob[100][5][5][5][5][5][5];
double x6prob[100][5][5][5][5][5][5];
double x7prob[100][5][5][5][5][5][5];
double x8prob[100][5][5][5][5][5][5];
double x9prob[100][5][5][5][5][5];
double x10prob[100][5][5][5][5][5][5];
double x11prob[100][5][5][5][5][5][5];
double x12prob[100][5][5][5][5][5][5];
double x13prob[100][5][5];
double x14prob[100][5][5][5][5][5][5];
double x15prob[100][5][5][5][5][5][5];
double x16prob[100][5][5][5][5][5][5];
double x17prob[100][5][5][5][5][5][5];
double x18prob[100][5][5][5][5][5][5];
double x19prob[100][5][5][5][5][5][5];
double x20prob[100][5][5][5][5][5][5];
double x21prob[100][5][5][5][5][5][5];
double x22prob[100][5][5][5][5][5];
double x23prob[100][5][5][5][5][5][5];
double x24prob[100][5][5][5][5][5];
double x25prob[100][5][5][5][5][5][5];
double x26prob[100][5][5][5][5][5];
double x27prob[100][5][5][5][5][5][5];
double x28prob[100][5][5][5][5][5][5];
double x29prob[100][5][5][5][5][5][5];
double x30prob[100][5][5];
double x31prob[100][5][5];
double x32prob[100][5][5];
double result[33][101];
double x[33][101][5];

int P[]={1,};
double weights[32];
int k[]={-1,2,2,2,2,0,0,0,0,0,0,2,2,0,0,2,0,2,0,0,0,0,0,2,0,0,0,2,2,2,0,0,0,2,2,0,0,2,2,2,2,0,0,2,2,0,0,0,0,};
int xbinNum[]={-1,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,};
double expdata[100][32];

void load0(char *path){

//loadExp
char buffer[MAX_LINE];
char *cp;
FILE *fin;

double tmp;
double sum[32];
int i,j,k;
for(k=0;k<32;k++){ sum[k]=0.0;}
fin = fopen("egfngf_cpu.csv","r");
fgets(buffer, sizeof(buffer), fin); //header
for(j=0;j<100;j++){
fgets(buffer, sizeof(buffer), fin);
cp = (char *)strtok(buffer, "\t, "); //time
for(k=0;k<32;k++){
cp = (char *)strtok(NULL, "\t, ");
tmp=strtod(cp, NULL);
expdata[j][k]=tmp;
sum[k]+=tmp;
}
}
fclose(fin);
for(k=0;k<32;k++){ weights[k]=100*1.0/sum[k];}

// load CPD
int idx;double p;
int vb,vi0,vi1,vi2,vi3,vi4,vi5,vi6,vi7,vi8,vi9,vi10;
int pi0,pi1,pi2,pi3,pi4,pi5,pi6,pi7,pi8,pi9,pi10;
for(tb=0;tb<100;tb++){

snprintf(buffer, sizeof(buffer), "%s/x0ctr_%d.txt",path, tb);
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
x1prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x1ctr_%d.txt",path, tb);
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
x2prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x2ctr_%d.txt",path, tb);
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

snprintf(buffer, sizeof(buffer), "%s/x3ctr_%d.txt",path, tb);
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
x4prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x4ctr_%d.txt",path, tb);
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
x5prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x5ctr_%d.txt",path, tb);
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
x6prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x6ctr_%d.txt",path, tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi4=idx%5;
idx=idx/5;
vi3=idx%5;
idx=idx/5;
vi2=idx%5;
idx=idx/5;
vi1=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
x7prob[tb][vi0][vi1][vi2][vi3][vi4][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x7ctr_%d.txt",path, tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi4=idx%5;
idx=idx/5;
vi3=idx%5;
idx=idx/5;
vi2=idx%5;
idx=idx/5;
vi1=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
x8prob[tb][vi0][vi1][vi2][vi3][vi4][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x8ctr_%d.txt",path, tb);
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
pi1=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x9prob[tb][pi0][pi1][vi0][vi1][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x9ctr_%d.txt",path, tb);
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
x10prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x10ctr_%d.txt",path, tb);
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
x11prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x11ctr_%d.txt",path, tb);
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
x12prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x12ctr_%d.txt",path, tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
x13prob[tb][vi0][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x13ctr_%d.txt",path, tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi3=idx%5;
idx=idx/5;
vi2=idx%5;
idx=idx/5;
vi1=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x14prob[tb][pi0][vi0][vi1][vi2][vi3][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x14ctr_%d.txt",path, tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi3=idx%5;
idx=idx/5;
vi2=idx%5;
idx=idx/5;
vi1=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x15prob[tb][pi0][vi0][vi1][vi2][vi3][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x15ctr_%d.txt",path, tb);
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
x16prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x16ctr_%d.txt",path, tb);
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
x17prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x17ctr_%d.txt",path, tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi3=idx%5;
idx=idx/5;
vi2=idx%5;
idx=idx/5;
vi1=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x18prob[tb][pi0][vi0][vi1][vi2][vi3][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x18ctr_%d.txt",path, tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi3=idx%5;
idx=idx/5;
vi2=idx%5;
idx=idx/5;
vi1=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x19prob[tb][pi0][vi0][vi1][vi2][vi3][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x19ctr_%d.txt",path, tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi3=idx%5;
idx=idx/5;
vi2=idx%5;
idx=idx/5;
vi1=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x20prob[tb][pi0][vi0][vi1][vi2][vi3][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x20ctr_%d.txt",path, tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi3=idx%5;
idx=idx/5;
vi2=idx%5;
idx=idx/5;
vi1=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x21prob[tb][pi0][vi0][vi1][vi2][vi3][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x21ctr_%d.txt",path, tb);
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
pi0=idx%5;
idx=idx/5;
x22prob[tb][pi0][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x22ctr_%d.txt",path, tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi3=idx%5;
idx=idx/5;
vi2=idx%5;
idx=idx/5;
vi1=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x23prob[tb][pi0][vi0][vi1][vi2][vi3][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x23ctr_%d.txt",path, tb);
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
pi1=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x24prob[tb][pi0][pi1][vi0][vi1][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x24ctr_%d.txt",path, tb);
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
x25prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x25ctr_%d.txt",path, tb);
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
pi1=idx%5;
idx=idx/5;
pi0=idx%5;
idx=idx/5;
x26prob[tb][pi0][pi1][vi0][vi1][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x26ctr_%d.txt",path, tb);
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
x27prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x27ctr_%d.txt",path, tb);
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
x28prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x28ctr_%d.txt",path, tb);
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
x29prob[tb][pi0][pi1][vi0][vi1][vi2][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x29ctr_%d.txt",path, tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
x30prob[tb][vi0][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x30ctr_%d.txt",path, tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
x31prob[tb][vi0][vb]=p;
}

fclose(fin);

snprintf(buffer, sizeof(buffer), "%s/x31ctr_%d.txt",path, tb);
fin = fopen(buffer, "r");
while (fgets(buffer, sizeof(buffer), fin)) {
cp = (char *)strtok(buffer, "\t, ");
idx=strtol(cp, NULL, 0);
cp = (char *)strtok(NULL, "\t, ");
p=strtod(cp, NULL);
vb=idx%5;
idx=idx/5;
vi0=idx%5;
idx=idx/5;
x32prob[tb][vi0][vb]=p;
}

fclose(fin);
}
}
void simulate(int kval[]){
int t,b,tb,i,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15;
double sum=0,margin=1,norm=0;

for(i=0;i<1;i++) k[P[i]]=kval[i];

//init
for(i=0;i<33;i++)
for(t=0;t<tps+1;t++)
for(b=0;b<5;b++)
x[i][t][b]=0;

for(i=0;i<33;i++)
for(t=0;t<tps+1;t++)
result[i][t]=0;


x[1][0][4]=1;
x[2][0][4]=1;
x[3][0][4]=1;
x[4][0][0]=1;
x[5][0][4]=1;
x[6][0][0]=1;
x[7][0][4]=1;
x[8][0][0]=1;
x[9][0][4]=1;
x[10][0][0]=1;
x[11][0][4]=1;
x[12][0][0]=1;
x[13][0][4]=1;
x[14][0][4]=1;
x[15][0][0]=1;
x[16][0][4]=1;
x[17][0][0]=1;
x[18][0][4]=1;
x[19][0][0]=1;
x[20][0][4]=1;
x[21][0][0]=1;
x[22][0][4]=1;
x[23][0][0]=1;
x[24][0][4]=1;
x[25][0][0]=1;
x[26][0][4]=1;
x[27][0][0]=1;
x[28][0][4]=1;
x[29][0][0]=1;
x[30][0][4]=1;
x[31][0][4]=1;
x[32][0][4]=1;

for(t=0;t<100;t++){
tb=t;
// normalization
for(i=1;i<33;i++){
norm=0; 
for(b=0;b<xbinNum[i];b++) norm+=x[i][t][b]; if(norm==0) norm=1;
for(b=0;b<xbinNum[i];b++) x[i][t][b]/=norm;
}


// normalization2
for(i=33;i<33;i++){
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
x[1][t+1][b]+=x1prob[tb][k[1]][k[2]][i0][i1][i2][b]*x[1][t][i0]*x[3][t][i1]*x[4][t][i2];
}
// x2
for(b=0;b<5;b++){
result[2][t]+=x[2][t][b]*M[2][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[2][t+1][b]+=x2prob[tb][k[3]][k[4]][i0][i1][i2][b]*x[2][t][i0]*x[5][t][i1]*x[6][t][i2];
}
// x3
for(b=0;b<5;b++){
result[3][t]+=x[3][t][b]*M[3][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[3][t+1][b]+=x3prob[tb][k[1]][k[2]][i0][i1][i2][b]*x[1][t][i0]*x[3][t][i1]*x[4][t][i2];
}
// x4
for(b=0;b<5;b++){
result[4][t]+=x[4][t][b]*M[4][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[4][t+1][b]+=x4prob[tb][k[1]][k[2]][i0][i1][i2][b]*x[1][t][i0]*x[3][t][i1]*x[4][t][i2];
}
// x5
for(b=0;b<5;b++){
result[5][t]+=x[5][t][b]*M[5][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[5][t+1][b]+=x5prob[tb][k[3]][k[4]][i0][i1][i2][b]*x[2][t][i0]*x[5][t][i1]*x[6][t][i2];
}
// x6
for(b=0;b<5;b++){
result[6][t]+=x[6][t][b]*M[6][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[6][t+1][b]+=x6prob[tb][k[3]][k[4]][i0][i1][i2][b]*x[2][t][i0]*x[5][t][i1]*x[6][t][i2];
}
// x7
for(b=0;b<5;b++){
result[7][t]+=x[7][t][b]*M[7][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
for(i3=0;i3<5;i3++)
for(i4=0;i4<5;i4++)
x[7][t+1][b]+=x7prob[tb][i0][i1][i2][i3][i4][b]*x[4][t][i0]*x[6][t][i1]*x[7][t][i2]*x[8][t][i3]*x[10][t][i4];
}
// x8
for(b=0;b<5;b++){
result[8][t]+=x[8][t][b]*M[8][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
for(i3=0;i3<5;i3++)
for(i4=0;i4<5;i4++)
x[8][t+1][b]+=x8prob[tb][i0][i1][i2][i3][i4][b]*x[4][t][i0]*x[6][t][i1]*x[7][t][i2]*x[8][t][i3]*x[10][t][i4];
}
// x9
for(b=0;b<5;b++){
result[9][t]+=x[9][t][b]*M[9][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
x[9][t+1][b]+=x9prob[tb][k[27]][k[28]][i0][i1][b]*x[9][t][i0]*x[21][t][i1];
}
// x10
for(b=0;b<5;b++){
result[10][t]+=x[10][t][b]*M[10][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[10][t+1][b]+=x10prob[tb][k[27]][k[28]][i0][i1][i2][b]*x[9][t][i0]*x[10][t][i1]*x[21][t][i2];
}
// x11
for(b=0;b<5;b++){
result[11][t]+=x[11][t][b]*M[11][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[11][t+1][b]+=x11prob[tb][k[11]][k[12]][i0][i1][i2][b]*x[8][t][i0]*x[11][t][i1]*x[12][t][i2];
}
// x12
for(b=0;b<5;b++){
result[12][t]+=x[12][t][b]*M[12][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[12][t+1][b]+=x12prob[tb][k[11]][k[12]][i0][i1][i2][b]*x[8][t][i0]*x[11][t][i1]*x[12][t][i2];
}
// x13
for(b=0;b<5;b++){
result[13][t]+=x[13][t][b]*M[13][b];
for(i0=0;i0<5;i0++)
x[13][t+1][b]+=x13prob[tb][i0][b]*x[13][t][i0];
}
// x14
for(b=0;b<5;b++){
result[14][t]+=x[14][t][b]*M[14][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
for(i3=0;i3<5;i3++)
x[14][t+1][b]+=x14prob[tb][k[15]][i0][i1][i2][i3][b]*x[12][t][i0]*x[14][t][i1]*x[15][t][i2]*x[25][t][i3];
}
// x15
for(b=0;b<5;b++){
result[15][t]+=x[15][t][b]*M[15][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
for(i3=0;i3<5;i3++)
x[15][t+1][b]+=x15prob[tb][k[15]][i0][i1][i2][i3][b]*x[12][t][i0]*x[14][t][i1]*x[15][t][i2]*x[25][t][i3];
}
// x16
for(b=0;b<5;b++){
result[16][t]+=x[16][t][b]*M[16][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[16][t+1][b]+=x16prob[tb][k[43]][k[44]][i0][i1][i2][b]*x[16][t][i0]*x[17][t][i1]*x[29][t][i2];
}
// x17
for(b=0;b<5;b++){
result[17][t]+=x[17][t][b]*M[17][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[17][t+1][b]+=x17prob[tb][k[43]][k[44]][i0][i1][i2][b]*x[16][t][i0]*x[17][t][i1]*x[29][t][i2];
}
// x18
for(b=0;b<5;b++){
result[18][t]+=x[18][t][b]*M[18][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
for(i3=0;i3<5;i3++)
x[18][t+1][b]+=x18prob[tb][k[17]][i0][i1][i2][i3][b]*x[15][t][i0]*x[17][t][i1]*x[18][t][i2]*x[19][t][i3];
}
// x19
for(b=0;b<5;b++){
result[19][t]+=x[19][t][b]*M[19][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
for(i3=0;i3<5;i3++)
x[19][t+1][b]+=x19prob[tb][k[17]][i0][i1][i2][i3][b]*x[15][t][i0]*x[17][t][i1]*x[18][t][i2]*x[19][t][i3];
}
// x20
for(b=0;b<5;b++){
result[20][t]+=x[20][t][b]*M[20][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
for(i3=0;i3<5;i3++)
x[20][t+1][b]+=x20prob[tb][k[23]][i0][i1][i2][i3][b]*x[19][t][i0]*x[20][t][i1]*x[21][t][i2]*x[31][t][i3];
}
// x21
for(b=0;b<5;b++){
result[21][t]+=x[21][t][b]*M[21][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
for(i3=0;i3<5;i3++)
x[21][t+1][b]+=x21prob[tb][k[23]][i0][i1][i2][i3][b]*x[19][t][i0]*x[20][t][i1]*x[21][t][i2]*x[31][t][i3];
}
// x22
for(b=0;b<5;b++){
result[22][t]+=x[22][t][b]*M[22][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[22][t+1][b]+=x22prob[tb][k[29]][i0][i1][i2][b]*x[4][t][i0]*x[12][t][i1]*x[22][t][i2];
}
// x23
for(b=0;b<5;b++){
result[23][t]+=x[23][t][b]*M[23][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
for(i3=0;i3<5;i3++)
x[23][t+1][b]+=x23prob[tb][k[29]][i0][i1][i2][i3][b]*x[4][t][i0]*x[12][t][i1]*x[22][t][i2]*x[23][t][i3];
}
// x24
for(b=0;b<5;b++){
result[24][t]+=x[24][t][b]*M[24][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
x[24][t+1][b]+=x24prob[tb][k[33]][k[34]][i0][i1][b]*x[23][t][i0]*x[24][t][i1];
}
// x25
for(b=0;b<5;b++){
result[25][t]+=x[25][t][b]*M[25][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[25][t+1][b]+=x25prob[tb][k[33]][k[34]][i0][i1][i2][b]*x[23][t][i0]*x[24][t][i1]*x[25][t][i2];
}
// x26
for(b=0;b<5;b++){
result[26][t]+=x[26][t][b]*M[26][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
x[26][t+1][b]+=x26prob[tb][k[37]][k[38]][i0][i1][b]*x[6][t][i0]*x[26][t][i1];
}
// x27
for(b=0;b<5;b++){
result[27][t]+=x[27][t][b]*M[27][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[27][t+1][b]+=x27prob[tb][k[37]][k[38]][i0][i1][i2][b]*x[6][t][i0]*x[26][t][i1]*x[27][t][i2];
}
// x28
for(b=0;b<5;b++){
result[28][t]+=x[28][t][b]*M[28][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[28][t+1][b]+=x28prob[tb][k[39]][k[40]][i0][i1][i2][b]*x[27][t][i0]*x[28][t][i1]*x[29][t][i2];
}
// x29
for(b=0;b<5;b++){
result[29][t]+=x[29][t][b]*M[29][b];
for(i0=0;i0<5;i0++)
for(i1=0;i1<5;i1++)
for(i2=0;i2<5;i2++)
x[29][t+1][b]+=x29prob[tb][k[39]][k[40]][i0][i1][i2][b]*x[27][t][i0]*x[28][t][i1]*x[29][t][i2];
}
// x30
for(b=0;b<5;b++){
result[30][t]+=x[30][t][b]*M[30][b];
for(i0=0;i0<5;i0++)
x[30][t+1][b]+=x30prob[tb][i0][b]*x[30][t][i0];
}
// x31
for(b=0;b<5;b++){
result[31][t]+=x[31][t][b]*M[31][b];
for(i0=0;i0<5;i0++)
x[31][t+1][b]+=x31prob[tb][i0][b]*x[31][t][i0];
}
// x32
for(b=0;b<5;b++){
result[32][t]+=x[32][t][b]*M[32][b];
for(i0=0;i0<5;i0++)
x[32][t+1][b]+=x32prob[tb][i0][b]*x[32][t][i0];
}
}
}

// output
void output(){
FILE *out;
char buffer[256];
snprintf(buffer, sizeof(buffer), "dummy.txt");
int idx=0,i,tb;

snprintf(buffer, sizeof(buffer), "brown_T_gpu.csv");
out=fopen(buffer, "w");
fprintf(out,"Time,");
for(i=1;i<33;i++)
fprintf(out,"x%d,",i);
fprintf(out,"\n");
for(tb=0;tb<tps;tb++){
fprintf(out,"%d,",tb);
for(i=1;i<33;i++)
fprintf(out,"%E,",result[i][tb]);
fprintf(out,"\n");
}
fclose(out);
snprintf(buffer, sizeof(buffer), "brown_T_gpu.csv");
out=fopen(buffer, "w");
fprintf(out,"Time,");
for(i=1;i<33;i++)
fprintf(out,"x%d,",i);
fprintf(out,"\n");
for(tb=0;tb<tps;tb++){
fprintf(out,"%f,",tb*1.0/dt);
for(i=1;i<33;i++)
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

snprintf(buffer, sizeof(buffer), "brown_M_gpu.csv");
out=fopen(buffer, "w");
for(t=0;t<tps;t++)
for(i=1;i<33;i++)
for(b=0;b<xbinNum[i];b++)
fprintf(out,"%d,%d,%d,%E\n",t,i,b,x[i][t][b]);
fclose(out);
}
double eval(int kval[]){
double err=0,value=0;int i,j;
simulate(kval);
for(i=0;i<100;i++){
for(j=0;j<32;j++){
err=result[j+1][i]-expdata[i][j];
value+=weights[j]*err*err;
}
}
return value;
}
int main(int argc, char *argv[]){
if(argc==1){
load0("../../sim/egfngf/tmp");
} else {
char path[128];
strcpy(path,argv[1]);
strcat(path,"/tmp/");
load0(path);
}
int i;
int k[]={2,};
printf("%E\n",eval(k));
outputMarginal();

}
