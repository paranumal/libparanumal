#include "boltzmann3D.h"

void boltzmannRampFunction3D(dfloat t, dfloat *ramp, dfloat *drampdt){
  // zero to SP at t=0
  *ramp = 0.5*(1+tanh(10.*(t-0.1)));
  *drampdt = 0.5*10*(1-pow(tanh(10.*(t-0.1)),2));  
  return;
}

