#define CONST(nr, min, max, bins) const unsigned int bins_const_ ## nr = bins;
#define REC(nr, cst) 
#define VAR(nr, min, max, bins, bininit) const unsigned int bins_var_ ## nr = bins;
#define AVAR(nr, min, max, bins) const unsigned int bins_var_ ## nr = bins;

#define AEQUATION(rest...)
#define EQUATION(rest...)

#include "model.cu"

#undef CONST
#undef REC
#undef VAR
#undef AVAR
#undef EQUATION
#undef AEQUATION

#define CONST(rest...)
#define REC(rest...)
#define VAR(rest...)
#define AVAR(rest...)

#define K(n) [bins_const_## n]
#define X(n) [bins_var_ ## n]
#define AX(n) [bins_var_ ## n]

#define AEQUATION(nr, formula, list) unsigned int x ## nr ## ctr[BLOCKS] list [bins_var_ ## nr];
#define EQUATION(nr, formula, list) unsigned int x ## nr ## ctr[BLOCKS] list [bins_var_ ## nr];

/*
        Generates the data structure for bins, used by both the host and GPU
*/

struct  xctr {
        #include "model.cu"
};

#undef K
#undef X
#undef AX
#undef CONST
#undef REC
#undef VAR
#undef AVAR
#undef EQUATION
#undef AEQUATION

