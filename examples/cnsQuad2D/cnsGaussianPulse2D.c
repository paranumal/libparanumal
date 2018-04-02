#include <math.h>
#include "mesh2D.h"

void cnsGaussianPulse2D(dfloat x, dfloat y, dfloat t,
			dfloat *u, dfloat *v, dfloat *p){

  // dudt = -dpdx
  // dvdt = -dpdy
  // dpdt = -dudx -dvdy

  // boundary conditions
  // dpdn = 0
  // u.n = 0

  *u = 0;
  *v = 0;
  *p = exp(-30*(x*x+y*y));

}
