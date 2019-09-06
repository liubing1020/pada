#define RK(result, integVar, formula, new_step_size, isOK) \
		real __attribute__((unused)) xx_delta, F1, F2, F3, F4, F5, F6, init_value; \
		real __attribute__((unused)) err1, err3, err4, err5, err6, realerr,current_h, xx, REL_TOL, ABS_TOL, tol; \
		current_h = ((real) (new_step_size)); REL_TOL = (real)(1E-6); ABS_TOL = (real)(1E-6); \
		xx = integVar; xx_delta = xx; init_value = (real)(formula); \
		realerr = ((real)0.0); \
		xx_delta = xx; F1 =  current_h * ((real)(init_value)); \
		xx_delta = xx + ((real)1.0/4) * F1; F2 =  current_h *((real)(formula)); \
		xx_delta = xx + ((real)(3.0/32)) * F1 + ((real)(9.0/32)) * F2 ; F3 =  current_h * ((real)(formula)); \
		xx_delta = xx + ((real)(1932.0/2197)) * F1 + ((real)(-7200.0/2197)) * F2 + ((real)(7296.0/2197)) * F3 ; \
		F4 = current_h * (real)(formula); \
		xx_delta = xx + ((real)(439.0/216)) * F1 + ((real)(-8.0)) * F2 + ((real)(3680.0/513)) * F3 + ((real)(-845.0/4104)) * F4 ; \
		F5 =  current_h * (real)(formula);\
		xx_delta = xx + ((real)(-8.0/27)) * F1 + ((real)2.0) * F2 + ((real)(-3544.0/2565)) * F3 + ((real)(1859.0/4104)) * F4 + ((real)(-11/40)) * F5; \
		F6 = current_h * (real)(formula); \
		err1 = ((real)(1.0/360)) * F1; \
		err3 = ((real)(-128.0/4275)) * F3; \
		err4 = ((real)(-2197.0/75240)) * F4; \
		err5 = ((real)(1.0/50)) * F5; \
		err6 = ((real)(2.0/55)) * F6; \
		realerr = (real)(err1 + err3 + err4 + err5 + err6) ; \
		tol = max(fabs((((real)REL_TOL) * init_value)), ABS_TOL);\
		if (fabs(((real)realerr)) > tol) {\
		new_step_size = ((real)0.84) * ((real)current_h * pow(fabs( tol / ((real)realerr)), ((real)0.25))); isOK = 0;} \
		else {isOK = 1; }\
		result = ((real)(25.0/216)) * F1 + ((real)(1408.0/2565)) * F3 + ((real)(2197.0/4104)) * F4 + ((real)(-1.0/5)) * F5 ; \
 

