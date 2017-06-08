#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh3D.h"

#define WADG 1

void acousticsSetup3D(mesh3D *mesh);

void acousticsRun3Dbbdg(mesh3D *mesh);
void acousticsOccaRun3Dbbdg(mesh3D *mesh);

void acousticsMRABUpdate3D(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat dt);  
void acousticsMRABUpdateTrace3D(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat dt); 
void acousticsMRABUpdate3D_wadg(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat dt); 
void acousticsMRABUpdateTrace3D_wadg(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat dt);  

void acousticsVolume3Dbbdg(mesh3D *mesh, iint lev);
void acousticsSurface3Dbbdg(mesh3D *mesh, iint lev, dfloat time);

void acousticsError3D(mesh3D *mesh, dfloat time);

void acousticsCavitySolution3D(dfloat x, dfloat y, dfloat z, dfloat time,
			       dfloat *u, dfloat *v, dfloat *w, dfloat *p);
void acousticsGaussianPulse3D(dfloat x, dfloat y, dfloat z, dfloat time,
			       dfloat *u, dfloat *v, dfloat *w, dfloat *p);

