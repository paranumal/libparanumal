#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh3D.h"

void acousticsRun3D(mesh3D *mesh);

void acousticsSetup3D(mesh3D *mesh);

void acousticsVolume3D(mesh3D *mesh);

void acousticsSurface3D(mesh3D *mesh, dfloat time);

void acousticsUpdate3D(mesh3D *mesh, dfloat rka, dfloat rkb);

void acousticsError3D(mesh3D *mesh, dfloat time);
