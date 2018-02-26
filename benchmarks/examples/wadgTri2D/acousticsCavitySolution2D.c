#include <math.h>
#include "mesh2D.h"

void acousticsCavitySolution2D(dfloat x, dfloat y, dfloat t,
			       dfloat *u, dfloat *v, dfloat *p){

  // dudt = -dpdx
  // dvdt = -dpdy
  // dpdt = -dudx -dvdy

  // boundary conditions
  // dpdn = 0
  // u.n = 0

  dfloat pi = M_PI;
  
  *u = sin(pi*t/sqrt(2.))*cos(pi*x/2.)*sin(pi*y/2.);
  *v = sin(pi*t/sqrt(2.))*sin(pi*x/2.)*cos(pi*y/2.);
  *p = (-sqrt(2.))*cos(pi*t/sqrt(2.))*sin(pi*x/2.)*sin(pi*y/2.);

}
