#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh3D.h"


void acousticsRun3Dbbdg(mesh3D *mesh);

void acousticsOccaRun3Dbbdg(mesh3D *mesh);

void acousticsSetup3D(mesh3D *mesh);

void acousticsVolume3Dbbdg(mesh3D *mesh);
void acousticsSurface3Dbbdg(mesh3D *mesh, dfloat time);
void acousticsUpdate3D(mesh3D *mesh, dfloat rka, dfloat rkb);

void acousticsError3D(mesh3D *mesh, dfloat time);

void acousticsCavitySolution3D(dfloat x, dfloat y, dfloat z, dfloat time,
			       dfloat *u, dfloat *v, dfloat *w, dfloat *p);

