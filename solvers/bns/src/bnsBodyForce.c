#include "bns.h"

void bnsBodyForce(dfloat t, dfloat *fx, dfloat *fy, dfloat *fz,
		  dfloat *intfx, dfloat *intfy, dfloat *intfz){

  *intfx = 0.5*(1+tanh(10.*(t-0.1)));
  *intfy = 0;
  *intfz = 0;

  *fx = 0.5*10*(1-pow(tanh(10.*(t-0.1)),2));
  *fy = 0;
  *fz = 0;
  
  return;
}

