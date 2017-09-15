#include <math.h>
#include "mesh3D.h"

void acousticsCavitySolution3D(dfloat x, dfloat y, dfloat z, dfloat time,
			       dfloat *u, dfloat *v, dfloat *w, dfloat *p){

  // dudt = -dpdx
  // dvdt = -dpdy
  // dwdt = -dpdz
  // dpdt = -dudx -dvdy -dwdz

  // boundary conditions
  // dpdn = 0
  // u.n = 0

  dfloat pi = M_PI;
  
  *u = sin(pi*time*sqrt(3.)/2.)*cos(pi*x/2.)*sin(pi*y/2.)*sin(pi*z/2.);
  *v = sin(pi*time*sqrt(3.)/2.)*sin(pi*x/2.)*cos(pi*y/2.)*sin(pi*z/2.);
  *w = sin(pi*time*sqrt(3.)/2.)*sin(pi*x/2.)*sin(pi*y/2.)*cos(pi*z/2.);
  *p = -sqrt(3.)*cos(pi*time*sqrt(3.)/2.)*sin(pi*x/2.)*sin(pi*y/2.)*sin(pi*z/2.);

}

