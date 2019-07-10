#define RK(result, integVar, formula) \
	real xx, __attribute__((unused)) xx_delta, F1, F2;\
	xx = integVar;\
	xx_delta = xx;\
	F1 = (real)(formula);\
	xx_delta = xx + ((real)(2.0/3)) * F1 * ((real)DT);\
	F2 = (real)(formula);\
	result = (((real)(1.0/4)) * F1 + ((real)(3.0/4)) * F2 ) * ((real) DT);
