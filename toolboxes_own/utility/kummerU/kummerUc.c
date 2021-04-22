#include <mex.h>
#include <matrix.h>
#include "arb.h"
#include "arb_hypgeom.h"

/* The computational routine 
void kummerUc(char *y, double *a_double, double *b_double, double *z_double, double *prec_double)
{	
	arb_t a, b, z, r;
	
	slong prec;
	prec = *prec_double;

	arb_init(a);
	arb_init(b);
	arb_init(z);
	arb_init(r);
	
	arb_set_d(a, *a_double);
	arb_set_d(b, *b_double);
	arb_set_d(z, *z_double);

	arb_hypgeom_u(r, a, b, z, prec);
	
	y = arb_get_str(r,1000,ARB_STR_NO_RADIUS);
	
	mexPrintf(y);
	
	arb_clear(a);
	arb_clear(b);
	arb_clear(z);
	arb_clear(r);
}
*/

/* The gateway function */
void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
	if ((nrhs<4)||(nrhs>4)) 
		mexErrMsgTxt("Must have four input arguments."); 
	
	if (nlhs>2)
		mexErrMsgTxt("Too many output arguments.");
		
	if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]))
		mexErrMsgTxt("All inputs must be double arrays.");
	
	double *a_double;
	double *b_double;
	double *z_double;
	double *prec_double;
	
	a_double = mxGetDoubles(prhs[0]);
	b_double = mxGetDoubles(prhs[1]);
	z_double = mxGetDoubles(prhs[2]);
	prec_double = mxGetDoubles(prhs[3]);
	
	/* Initialize arb_t. arb_t represents a ball over the real numbers, that is, an interval [m±r]≡[m−r,m+r] where the midpoint m and the radius r are (extended) real numbers and r is nonnegative (possibly infinite) */
	arb_t a, b, z, r;
	
	/* Set precision parameter*/
	slong prec;
	prec = *prec_double;

	arb_init(a);
	arb_init(b);
	arb_init(z);
	arb_init(r);
	
	arb_set_d(a, *a_double);
	arb_set_d(b, *b_double);
	arb_set_d(z, *z_double);

	arb_hypgeom_u(r, a, b, z, prec);
	
	char *y;
	/* This converts the ball r, to a string which is just the midpoint of the ball */
	y = arb_get_str(r,1000,ARB_STR_NO_RADIUS);
	
	/*mexPrintf(y);*/
	
	/*char *y;*/
    /* allocate memory for output string 
    y=mxCalloc(1000, sizeof(char));	*/
	
	/* call computational routine */
	/*kummerUc(y,a,b,z,prec);*/
	
	/* output the string y */
	plhs[0] = mxCreateString(y);
	
	/*
	mxFree(a_double);
	mxFree(b_double);
	mxFree(z_double);
	mxFree(prec_double);
	*/
	
	arb_clear(a);
	arb_clear(b);
	arb_clear(z);
	arb_clear(r);
	
	flint_cleanup();
	
	return;
}