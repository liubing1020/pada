		#define RK(result, integVar, formula, new_step_size, isOK) \
		real __attribute__((unused)) xx_delta, F1, F2, F3, F4, F5, F6, init_value; \
		real __attribute__((unused)) err1, err3, err4, err5, err6, realerr,current_h, xx, REL_TOL, ABS_TOL, tol; \
		current_h = ((real) (new_step_size)); REL_TOL = (real)(1E-6); ABS_TOL = (real)(1E-6); \
		xx = integVar; xx_delta = xx; init_value = (real)(formula); \
		realerr = ((real)0.0); \
		xx_delta = xx; F1 =  current_h * ((real)(init_value)); \
		xx_delta = xx + ((real)1.0/5) * F1; F2 =  current_h *((real)(formula)); \
		xx_delta = xx + ((real)(3.0/40)) * F1 + ((real)(9.0/40)) * F2 ; F3 =  current_h * ((real)(formula)); \
		xx_delta = xx + ((real)(3.0/10)) * F1 + ((real)(-9.0/10)) * F2 + ((real)(6.0/5)) * F3 ; \
		F4 = current_h * (real)(formula); \
		xx_delta = xx + ((real)(-11.0/54)) * F1 + ((real)(5.0/2)) * F2 + ((real)(-70.0/27)) * F3 + ((real)(35.0/27)) * F4 ; \
		F5 =  current_h * (real)(formula);\
		xx_delta = xx + ((real)(1631.0/55296)) * F1 + ((real)(175/512)) * F2 + ((real)(575.0/13824)) * F3 + ((real)(44275.0/110592)) * F4 + ((real)(253.0/4096)) * F5 ; \
		F6 = current_h * (real)(formula); \
		err1 = ((real)(-277.0/64512)) * F1; \
		err3 = ((real)(0.018668586)) * F3; \
		err4 = ((real)(-0.034155026)) * F4; \
		err5 = ((real)(-277.0/14336)) * F5; \
		err6 = ((real)(277.0/7084)) * F6; \
		realerr = (real)(err1 + err3 + err4 + err5 + err6) ; \
		tol = max(fabs((((real) REL_TOL) * init_value)), ABS_TOL);\
		if (fabs(((real)realerr)) > tol) {\
		new_step_size = ((real)0.84) * ((real)current_h * pow(fabs( tol / ((real)realerr)), ((real)0.25))); \
		isOK = 0;} \
		else {isOK = 1; }\
		result = ((real)(37.0/378)) * F1 + ((real)(250.0/621)) * F3 + ((real)(125.0/594)) * F4 + ((real)(512.0/1771)) * F6 ; \
 

