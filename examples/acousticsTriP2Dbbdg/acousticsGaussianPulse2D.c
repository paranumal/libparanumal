#include <math.h>
#include "mesh2D.h"

void acousticsGaussianPulse2D(dfloat x, dfloat y, dfloat t,
			      dfloat *u, dfloat *v, dfloat *p){

  *u = 0;
  *v = 0;
  *p = exp(-30*(x*x+y*y));

}
