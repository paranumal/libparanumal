#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh3D.h"

#define WADG 1
#define ASYNC 1

void acousticsSetup3D(mesh3D *mesh);

void acousticsRun3Dbbdg(mesh3D *mesh);
void acousticsOccaRun3Dbbdg(mesh3D *mesh);

void acousticsMRABUpdate3D(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat t, dfloat dt);  
void acousticsMRABUpdateTrace3D(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat t, dfloat dt); 
void acousticsMRABUpdate3D_wadg(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat t, dfloat dt); 
void acousticsMRABUpdateTrace3D_wadg(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat t, dfloat dt);  

void acousticsVolume3Dbbdg(mesh3D *mesh, iint lev);
void acousticsSurface3Dbbdg(mesh3D *mesh, iint lev, dfloat time);

void acousticsPmlSetup3D(mesh3D *mesh);

void acousticsPmlVolume3Dbbdg(mesh3D *mesh, iint lev);
void acousticsPmlSurface3Dbbdg(mesh3D *mesh, iint lev, dfloat t);
void acousticsMRABpmlUpdate3D(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat dt);
void acousticsMRABpmlUpdateTrace3D(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat dt);
void acousticsMRABpmlUpdate3D_wadg(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat dt);
void acousticsMRABpmlUpdateTrace3D_wadg(mesh3D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat dt);

void acousticsSourceSetup3D(mesh3D *mesh, occa::kernelInfo &kernelInfo);

void acousticsError3D(mesh3D *mesh, dfloat time);

void acousticsCavitySolution3D(dfloat x, dfloat y, dfloat z, dfloat time,
			       dfloat *u, dfloat *v, dfloat *w, dfloat *p);
void acousticsGaussianPulse3D(dfloat x, dfloat y, dfloat z, dfloat time,
			       dfloat *u, dfloat *v, dfloat *w, dfloat *p);

