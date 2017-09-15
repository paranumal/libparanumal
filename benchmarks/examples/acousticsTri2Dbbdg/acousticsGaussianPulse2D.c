#include <math.h>
#include "mesh2D.h"

void acousticsGaussianPulse2D(dfloat x, dfloat y, dfloat t,
			      dfloat *u, dfloat *v, dfloat *p){

  *u = 0;
  *v = 0;
  dfloat x0 = x - 0.;
  dfloat y0 = y - 0.2;

  *p = exp(-30*(x0*x0+y0*y0));

}
