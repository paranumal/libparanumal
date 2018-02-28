#include <math.h>
#include "mesh3D.h"

void acousticsGaussianPulse3D(dfloat x, dfloat y, dfloat z, dfloat time,
			       dfloat *u, dfloat *v, dfloat *w, dfloat *p){

  *u = 0.;
  *v = 0.;
  *w = 0.;
  dfloat x0 = x - 0.;
  dfloat y0 = y - 0.;
  dfloat z0 = z - 0.2;

  *p = exp(-30*(x0*x0+y0*y0+z0*z0));

}