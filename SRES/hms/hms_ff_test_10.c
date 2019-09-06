#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

void load9(char*);
void simulate(int kval[]);
void output(void);
void outputMarginal(void);
double eval(int kval[]);

int tps=40,tb,MAX_LINE=1024;
float dt=1;
double M[][16]={{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,},
		{179.0077753229652,225.35743726836873,283.7082045225821,357.1674682202972,449.647202031126,566.0722890195802,712.6427895550448,897.164117286746,1129.462705813709,1421.9093020347073,1790.077753598196,2253.5743730522313,2837.082045594365,3571.6746825715163,4496.472020679805,5660.722890564346,},
		{179.0077753229652,225.35743726836873,283.7082045225821,357.1674682202972,449.647202031126,566.0722890195802,712.6427895550448,897.164117286746,1129.462705813709,1421.9093020347073,1790.077753598196,2253.5743730522313,2837.082045594365,3571.6746825715163,4496.472020679805,5660.722890564346,},
		{3.810988127367627E-5,1.7908962108401907E-4,5.596351460205221E-4,0.0015868385489287724,0.004359561031915262,0.011843951380065793,0.032046512836079406,0.08657914133996075,0.2337786744922155,0.6311133180750679,1.7036358699954277,4.5986882973982,12.41328253991852,33.50716069331447,90.44571459392668,244.13954451755035,},
};

double x1prob[40][5][5][5][5][5][16][16][16];
double x2prob[40][5][5][5][5][16][16][16];
double x3prob[40][5][5][5][5][5][16][16][16];
double result[4][41];
double x[4][41][16];

int P[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,};
double weights[3];
int k[]={-1,2,2,2,2,2,2,2,2,2,2,2,2,4,4,};
int xbinNum[]={-1,16,16,16,};
double expdata[40][3];

void load9(char *path){

  //loadExp
  char buffer[MAX_LINE];
  char *cp;
  FILE *fin;

  double tmp;
  double sum[3];
  int i,j,k;
  for(k=0;k<2;k++){ sum[k]=0.0;}
  fin = fopen("./40_points.txt","r");
  fgets(buffer, sizeof(buffer), fin); //header
  for(j=0;j<40;j++){
    fgets(buffer, sizeof(buffer), fin);
    cp = (char *)strtok(buffer, "\t, "); //time
    for(k=0;k<2;k++){
      cp = (char *)strtok(NULL, "\t, ");
      tmp=strtod(cp, NULL);
      expdata[j][k]=tmp;
      sum[k]+=tmp;
    }
  }
  fclose(fin);
  for(k=0;k<2;k++){ weights[k]=40*1.0/sum[k];}

  // load CPD
  int idx;double p;
  int vb,vi0,vi1,vi2,vi3,vi4,vi5,vi6,vi7,vi8,vi9,vi10;
  int pi0,pi1,pi2,pi3,pi4,pi5,pi6,pi7,pi8,pi9,pi10;
  for(tb=0;tb<40;tb++){ //<--

    snprintf(buffer, sizeof(buffer), "%s/x0ctr_%d.txt",path, tb);
    fin = fopen(buffer, "r");
    while (fgets(buffer, sizeof(buffer), fin)) {
      cp = (char *)strtok(buffer, "\t, ");
      idx=strtol(cp, NULL, 0);
      cp = (char *)strtok(NULL, "\t, ");
      p=strtod(cp, NULL);
      vb=idx%16;
      idx=idx/16;
      vi1=idx%16;
      idx=idx/16;
      vi0=idx%16;
      idx=idx/16;
      pi4=idx%5;
      idx=idx/5;
      pi3=idx%5;
      idx=idx/5;
      pi2=idx%5;
      idx=idx/5;
      pi1=idx%5;
      idx=idx/5;
      pi0=idx%5;
      idx=idx/5;
      x1prob[tb][pi0][pi1][pi2][pi3][pi4][vi0][vi1][vb]=p; 
	  //printf("%d,%d,%d,%d,%d,%d,%d,%d\n",pi0,pi1,pi2,pi3,pi4,vi0,vi1,vb);
    }

	//printf("\n");
    fclose(fin);

    snprintf(buffer, sizeof(buffer), "%s/x1ctr_%d.txt",path, tb);
    fin = fopen(buffer, "r");
    while (fgets(buffer, sizeof(buffer), fin)) {
      cp = (char *)strtok(buffer, "\t, ");
      idx=strtol(cp, NULL, 0);
      cp = (char *)strtok(NULL, "\t, ");
      p=strtod(cp, NULL);
      vb=idx%16;
      idx=idx/16;
      vi1=idx%16;
      idx=idx/16;
      vi0=idx%16;
      idx=idx/16;
      pi3=idx%5;
      idx=idx/5;
      pi2=idx%5;
      idx=idx/5;
      pi1=idx%5;
      idx=idx/5;
      pi0=idx%5;
      idx=idx/5;
      x2prob[tb][pi0][pi1][pi2][pi3][vi0][vi1][vb]=p;
	  //printf("%d,%d,%d,%d,%d,%d,%d\n",pi0,pi1,pi2,pi3,vi0,vi1,vb);
    }

		//printf("\n");
    fclose(fin);

    snprintf(buffer, sizeof(buffer), "%s/x2ctr_%d.txt",path, tb);
    fin = fopen(buffer, "r");
    while (fgets(buffer, sizeof(buffer), fin)) {
      cp = (char *)strtok(buffer, "\t, ");
      idx=strtol(cp, NULL, 0);
      cp = (char *)strtok(NULL, "\t, ");
      p=strtod(cp, NULL);
      vb=idx%16;
      idx=idx/16;
      vi1=idx%16;
      idx=idx/16;
      vi0=idx%16;
      idx=idx/16;
      pi4=idx%5;
      idx=idx/5;
      pi3=idx%5;
      idx=idx/5;
      pi2=idx%5;
      idx=idx/5;
      pi1=idx%5;
      idx=idx/5;
      pi0=idx%5;
      idx=idx/5;
      x3prob[tb][pi0][pi1][pi2][pi3][pi4][vi0][vi1][vb]=p;
	  //printf("%d,%d,%d,%d,%d,%d,%d,%d\n",pi0,pi1,pi2,pi3,pi4,vi0,vi1,vb);
    }

    fclose(fin);
  }
}
void simulate(int kval[]){
  char buffer[MAX_LINE];
  char *cp;
  FILE *fin;
  int t,b,tb,j,i,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15;
  double sum=0,margin=1,norm=0;

  for(i=0;i<14;i++) k[P[i]]=kval[i];

  //init
  for(i=0;i<4;i++)
    for(t=0;t<tps+1;t++)
      for(b=0;b<16;b++)
	x[i][t][b]=0;

  for(i=0;i<4;i++)
    for(t=0;t<tps+1;t++)
      result[i][t]=0;


  x[1][0][9]=1;
  x[2][0][12]=1;

  for(t=0;t<40;t++){
    tb=t;
    // normalization
    for(i=1;i<3;i++){
      norm=0; 
      for(b=0;b<xbinNum[i];b++) norm+=x[i][t][b]; if(norm==0) norm=1;
      for(b=0;b<xbinNum[i];b++) x[i][t][b]/=norm;
    }

    // x3
    for(b=0;b<16;b++){
      for(i0=0;i0<16;i0++)
	for(i1=0;i1<16;i1++)
	  x[3][t][b]+=x3prob[tb][k[4]][k[5]][k[6]][k[12]][k[14]][i0][i1][b]*x[1][t][i0]*x[2][t][i1];
      result[3][t]+=x[3][t][b]*M[3][b];
    }

    // normalization2
    for(i=3;i<4;i++){
      norm=0; 
      for(b=0;b<xbinNum[i];b++) norm+=x[i][t][b];if(norm==0) norm=1;
      for(b=0;b<xbinNum[i];b++) x[i][t][b]/=norm;
    }

    // x1
    for(b=0;b<16;b++){
      result[1][t]+=x[1][t][b]*M[1][b];
      for(i0=0;i0<16;i0++)
	for(i1=0;i1<16;i1++)
	  x[1][t+1][b]+=x1prob[tb][k[1]][k[2]][k[3]][k[11]][k[13]][i0][i1][b]*x[1][t][i0]*x[3][t][i1];
    }
    // x2
    for(b=0;b<16;b++){
      result[2][t]+=x[2][t][b]*M[2][b];
      for(i0=0;i0<16;i0++)
	for(i1=0;i1<16;i1++)
	  x[2][t+1][b]+=x2prob[tb][k[7]][k[8]][k[9]][k[10]][i0][i1][b]*x[1][t][i0]*x[2][t][i1];
    }
  }
}

// output
void output(){
  FILE *out;
  char buffer[256];
  snprintf(buffer, sizeof(buffer), "dummy.txt");
  int idx=0,i,tb;

  snprintf(buffer, sizeof(buffer), "hms_T_gpu.csv");
  out=fopen(buffer, "w");
  fprintf(out,"Time,");
  for(i=1;i<4;i++)
    fprintf(out,"x%d,",i);
  fprintf(out,"\n");
  for(tb=0;tb<tps;tb++){
    fprintf(out,"%d,",tb);
    for(i=1;i<4;i++)
      fprintf(out,"%E,",result[i][tb]);
    fprintf(out,"\n");
  }
  fclose(out);
  snprintf(buffer, sizeof(buffer), "hms_T_gpu.csv");
  out=fopen(buffer, "w");
  fprintf(out,"Time,");
  for(i=1;i<4;i++)
    fprintf(out,"x%d,",i);
  fprintf(out,"\n");
  for(tb=0;tb<tps;tb++){
    fprintf(out,"%f,",tb*1.0/dt);
    for(i=1;i<4;i++)
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

  snprintf(buffer, sizeof(buffer), "hms_M_gpu.csv");
  out=fopen(buffer, "w");
  for(t=0;t<tps;t++)
    for(i=1;i<4;i++)
      for(b=0;b<xbinNum[i];b++)
	fprintf(out,"%d,%d,%d,%E\n",t,i,b,x[i][t][b]);
  fclose(out);
}
double eval(int kval[]){
  double err=0,value=0;int i,j;
  simulate(kval);
  for(i=0;i<40;i++){
    for(j=0;j<2;j++){
      err=result[j+1][i]-expdata[i][j];
      value+=weights[j]*err*err;
    }
  }
  return value;
}


int main(int argc, char *argv[]){
  if(argc==1){
    load9("/hpctmp/dcsliubi/hms/6M/10/tmp/");
  } else {
    char path[128];
    strcpy(path,argv[1]);
    strcat(path,"/tmp/");
    load9(path);
  }
  struct timeval start = {0,0};
  struct timeval end = {0,0};
  float seconds;
  srand ( time(NULL) );
  int i,j;
  //int k[]={2,2,2,2,2,2,2,2,2,2,2,2,4,4,};
  //int k[]={4,0,3,2,4,0,0,0,4,3,2,1,2,4,};
  //int k[]={4,2,3,2,4,1,0,0,4,0,0,4,1,4,}; // results dec 13
  //int k[]={1,2,2,4,0,1,0,0,4,0,0,4,4,2}; // results dec 15
  //int k[]={4,3,2,0,4,0,0,0,4,0,0,4,2,4}; // results dec 16
  int k[]={1,4,4,3,3,0,0,0,4,0,0,4,3,4,};
  gettimeofday(&start, NULL);
  /*
  for(i=0;i<1000;i++){ 
	for(j=0;j<14;j++){
		k[j]=rand() % 5;
		printf("%d",k[j]);
	}
    printf("          %E\n",eval(k));
  }
*/
printf("          %E\n",eval(k));
output();

  gettimeofday(&end, NULL);
  seconds = (end.tv_sec-start.tv_sec) + (float)(end.tv_usec-start.tv_usec) / 1000000;
  printf("elapsed time\n> %f seconds\n", seconds);
  printf("%E\n",eval(k));


}

