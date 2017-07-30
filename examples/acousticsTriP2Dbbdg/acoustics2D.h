
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"

#define WADG 1

void acousticsSetup2D(mesh2D *mesh);

void acousticsRun2Dbbdg(mesh2D *mesh);
void acousticsOccaRun2Dbbdg(mesh2D *mesh);

void acousticsMRABUpdate2D(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat t, dfloat dt);  
void acousticsMRABUpdateTrace2D(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat t, dfloat dt); 
void acousticsMRABUpdate2D_wadg(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat t, dfloat dt); 
void acousticsMRABUpdateTrace2D_wadg(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat t, dfloat dt);  

void acousticsVolume2Dbbdg(mesh2D *mesh, iint lev);
void acousticsSurface2Dbbdg(mesh2D *mesh, iint lev, dfloat t);

void acousticsPmlSetup2D(mesh2D *mesh);

void acousticsPmlVolume2Dbbdg(mesh2D *mesh, iint lev);
void acousticsPmlSurface2Dbbdg(mesh2D *mesh, iint lev, dfloat t);
void acousticsMRABpmlUpdate2D(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat dt);
void acousticsMRABpmlUpdateTrace2D(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat dt);
void acousticsMRABpmlUpdate2D_wadg(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat dt);
void acousticsMRABpmlUpdateTrace2D_wadg(mesh2D *mesh, dfloat ab1, dfloat ab2, dfloat ab3, iint lev, dfloat dt);

void acousticsSourceSetup2D(mesh2D *mesh, occa::kernelInfo &kernelInfo);


void acousticsError2D(mesh2D *mesh, dfloat time);

void acousticsGaussianPulse2D(dfloat x, dfloat y, dfloat t, dfloat *u, dfloat *v, dfloat *p);
