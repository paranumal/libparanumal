#include <math.h>
#include "mesh3D.h"

void cnsGaussianPulse3D(dfloat x, dfloat y, dfloat z, dfloat t,
			dfloat *r, dfloat *u, dfloat *v, dfloat *w){

  *r = 1 + exp(-3*(x*x+y*y+z*z));
  *u = 0;
  *v = 0;
  *w = 0;

}
