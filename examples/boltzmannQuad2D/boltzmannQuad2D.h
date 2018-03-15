
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "mesh2D.h"


void boltzmannErrorQuad2D(mesh2D *mesh, dfloat time);

void boltzmannComputeVorticityQuad2D(mesh2D *mesh, dfloat *q, int outfld, int Nfields);

void boltzmannSplitPmlRunQuad2D(mesh2D *mesh);

void boltzmannSplitPmlSetupQuad2D(mesh2D *mesh);

dfloat boltzmannRampFunction2D(dfloat t);
