/*
(c) Andrei Hagiescu 2010
The file writer code, it is run on the host
*/

#define CONST(rest...)
#define REC(rest...)
#define VAR(rest...)
#define AVAR(rest...)

#define K(n) for (int k ## n = 0; k ## n < bins_const_ ## n; k ## n ++)

#define AX X
#define X(n) for (int x ## n = 0; x ## n < bins_var_ ## n; x ## n ++)

#define AEQUATION EQUATION
#define EQUATION(nr, formula, list) { \
	char buffer[256]; \
	FILE * out; \
	snprintf(buffer, sizeof(buffer), "tmp/x" #nr "ctr_%d.txt", block); \
	out = fopen(buffer, "w"); \
	int *structptr = (int *)(&bins_host-> x ## nr ## ctr[block]); \
	int idx = 0; \
	list {\
		real sum = 0; \
		for (int vi = 0; vi < bins_var_ ## nr; vi ++) sum += *(structptr + vi); \
		for (int vi = 0; vi < bins_var_ ## nr; vi ++) { \
			checkup_counter += *(structptr); \
			real elem = *(structptr++); \
			if (elem > 0 && sum > 0) { \
				elem /= sum; \
				fprintf(out, "%d %10.20f\n", idx, (double)elem); \
			} \
			idx++; \
		} \
	} \
	fclose(out); \
} \


for (int block = 0; block < BLOCKS; block++) {
		#include "model.cu"
}

#undef CONST
#undef REC
#undef VAR
#undef AVAR
#undef K
#undef X 
#undef AX
#undef EQUATION
#undef AEQUATION

