#include "arb.h"
#include "arb_hypgeom.h"
#include <stdio.h>

int main(int argc, char ** argv)
{	
	double a_double;
	double b_double;
	double z_double;
	double prec_double;
	
	a_double = atof(argv[1]);
	b_double = atof(argv[2]);
	z_double = atof(argv[3]);
	prec_double = atof(argv[4]);
	
	arb_t a, b, z, r;
	
	slong prec;
	prec = prec_double;

	arb_init(a);
	arb_init(b);
	arb_init(z);
	arb_init(r);
	
	arb_set_d(a, a_double);
	arb_set_d(b, b_double);
	arb_set_d(z, z_double);

	arb_hypgeom_u(r, a, b, z, prec);
	arb_printd(r,100);
	flint_printf("\n");
	
	char *y;
	
	y = arb_get_str(r,1000,ARB_STR_NO_RADIUS);
	
	arb_clear(a);
	arb_clear(b);
	arb_clear(z);
	arb_clear(r);
	
	printf(y);
	
	return 0;
}
