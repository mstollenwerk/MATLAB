#include <mex.h>
#include <matrix.h>
#include "arb.h"
#include "arb_hypgeom.h"

/* The gateway function */
void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
	if ((nrhs<3)||(nrhs>3)) 
		mexErrMsgTxt("Must have three input arguments."); 
	
	if (nlhs>2)
		mexErrMsgTxt("Too many output arguments.");
		
	if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]))
		mexErrMsgTxt("All inputs must be double arrays.");
	
	double *a_double;
	double *z_double;
	double *prec_double;
	
	a_double = mxGetDoubles(prhs[0]);
	z_double = mxGetDoubles(prhs[1]);
	prec_double = mxGetDoubles(prhs[2]);
	
	/* Initialize arb_t. arb_t represents a ball over the real numbers, that is, an interval [m±r]≡[m−r,m+r] where the midpoint m and the radius r are (extended) real numbers and r is nonnegative (possibly infinite) */
	arb_t a, z, r;
	
	/* Set precision parameter*/
	slong prec;
	prec = *prec_double;

	arb_init(a);
	arb_init(z);
	arb_init(r);
	
	arb_set_d(a, *a_double);
	arb_set_d(z, *z_double);

	arb_hypgeom_bessel_k(r, a, z, prec);
	
	char *y;
	/* This converts the ball r, to a string which is just the midpoint of the ball */
	y = arb_get_str(r,1000,ARB_STR_NO_RADIUS);

	plhs[0] = mxCreateString(y);
	
	arb_clear(a);
	arb_clear(z);
	arb_clear(r);
	
	flint_cleanup();
	
	return;
}