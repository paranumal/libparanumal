#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh3D.h"

void acousticsRunHex3D(mesh3D *mesh);

void acousticsOccaRunHex3D(mesh3D *mesh);

void acousticsSetupHex3D(mesh3D *mesh);

void acousticsVolumeHex3D(mesh3D *mesh);

void acousticsSurfaceHex3D(mesh3D *mesh, dfloat time);

void acousticsUpdateHex3D(mesh3D *mesh, dfloat rka, dfloat rkb);

void acousticsErrorHex3D(mesh3D *mesh, dfloat time);


