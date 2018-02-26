#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

#include "mesh2D.h"

void acousticsRunQuad2D(mesh2D *mesh);

void acousticsSetupQuad2D(mesh2D *mesh);

void acousticsVolumeQuad2D(mesh2D *mesh);

void acousticsSurfaceQuad2D(mesh2D *mesh);

void acousticsUpdate2D(mesh2D *mesh, dfloat rka, dfloat rkb);

void acousticsError2D(mesh2D *mesh, dfloat time);
