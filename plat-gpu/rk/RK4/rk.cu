#define RK(result, integVar, formula) \
	real xx, __attribute__((unused)) xx_delta, halfF1, halfF2, F3, F4;\
	xx = integVar;\
	xx_delta = xx;\
	halfF1 = ((real)HALFDT) * (formula);\
	xx_delta = xx + halfF1;\
	halfF2 = ((real)HALFDT) * (formula);\
	xx_delta = xx + halfF2;\
	F3 = ((real)DT) * (formula);\
	xx_delta = xx + F3;\
	F4 = ((real)DT) * (formula);\
	result = (((real)2.0) * (halfF1 + F3) + (real)4.0 * halfF2 + F4) * ((real)1.0/(real)6.0);

