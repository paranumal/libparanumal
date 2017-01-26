#include <math.h>
#include "mesh2D.h"

void boltzmannGaussianPulse2D(dfloat x, dfloat y, dfloat t,
			      dfloat *q1, dfloat *q2, dfloat *q3,
			      dfloat *q4, dfloat *q5, dfloat *q6){

  *q1 = 1. + .1*exp(-30*(x*x+y*y));

}
