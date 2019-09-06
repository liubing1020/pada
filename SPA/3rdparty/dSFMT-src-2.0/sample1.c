#include <stdio.h>
#include <stdlib.h>
#include "dSFMT.c"

int main(int argc, char* argv[]) {
    int i, cnt, seed;
    double x, y, pi;
    const int NUM = 100000000;
    dsfmt_t dsfmt;

    if (argc >= 2) {
	seed = strtol(argv[1], NULL, 10);
    } else {
	seed = 12345;
    }
    cnt = 0;
    dsfmt_init_gen_rand(&dsfmt,seed);
    x=0;
    for (i = 0; i < NUM; i++) {
      //printf("%f\n",genrand_real3());
      x = x+dsfmt_genrand_close_open(&dsfmt);
	//y = genrand_res53();
	//if (x * x + y * y < 1.0) {
	//  cnt++;
	//}
    }
    //pi = (double)cnt / NUM * 4;
    //printf("%lf\n", pi);
    printf("%f\n",x);
    return 0;
}
