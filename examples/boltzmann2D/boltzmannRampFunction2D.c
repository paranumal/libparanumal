#include "boltzmann2D.h"

dfloat boltzmannRampFunction2D(dfloat t){

  // zero to SP at t=0
  dfloat ramp = 0.5*(1+tanh(10.*(t-0.1)));

  return ramp;
}

