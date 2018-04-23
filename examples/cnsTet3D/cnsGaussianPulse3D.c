#include <math.h>
#include "mesh2D.h"

void cnsGaussianPulse2D(dfloat x, dfloat y, dfloat z, dfloat t,
			dfloat *r, dfloat *u, dfloat *v, dfloat *w){

  *r = 1 + exp(-30*(x*x+y*y+z*z));
  *u = 0;
  *v = 0;
  *w = 0;

}
