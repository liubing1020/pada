#define RK(result, integVar, formula) \
	real __attribute__((unused)) xx_delta, F1, F2, F3, xx; \
	xx = integVar; \
	xx_delta = xx; F1 = (real)(formula); \
	xx_delta = xx + ((real)(1.0/2)) * F1 * ((real)DT); F2 = (real)(formula); \
	xx_delta = xx + ((real)(3.0/4)) * F2 * ((real)DT); F3 = (real)(formula); \
	result  = (((real)(2.0/9)) * F1 + ((real)(1.0/3)) * F2 + ((real)(4.0/9)) * F3 ) * ((real) DT);

