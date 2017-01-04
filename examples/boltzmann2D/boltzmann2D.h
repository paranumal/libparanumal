
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "mesh2D.h"


void boltzmannError2D(mesh2D *mesh, dfloat time);

void boltzmannComputeVorticity2D(mesh2D *mesh, dfloat *q, iint outfld, iint Nfields);

//dfloat boltzmannRampFunction2D(dfloat t);
void boltzmannRampFunction2D(dfloat t, dfloat *ramp, dfloat *drampdt);
