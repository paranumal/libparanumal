#include "bns.h"

void bnsRampFunction(dfloat t, dfloat *ramp, dfloat *drampdt){

  // zero to SP at t=0
 #if 0	
  *ramp = 0.5*(1+tanh(10.*(t-0.1)));
  *drampdt = 0.5*10*(1-pow(tanh(10.*(t-0.1)),2));

 #else
   *ramp    = 1.0;
   *drampdt = 0.0;
 #endif
  
  return;
}

