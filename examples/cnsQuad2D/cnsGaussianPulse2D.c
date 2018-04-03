#include <math.h>
#include "mesh2D.h"

void cnsGaussianPulse2D(dfloat x, dfloat y, dfloat t,
			dfloat *r, dfloat *u, dfloat *v){

  *r = 1 + exp(-30*(x*x+y*y));
  *u = 0;
  *v = 0;

}
