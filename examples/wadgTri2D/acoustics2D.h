#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh2D.h"

void acousticsSetup2D(mesh2D *mesh);

void acousticsPmlSetup2D(mesh2D *mesh,
			 dfloat xmin, dfloat xmax, // bounding box for non-pml sub-domain
			 dfloat ymin, dfloat ymax,
			 dfloat xsigma, dfloat ysigma);

void acousticsSplitPmlSetup2D(mesh2D *mesh);

void acousticsPml2D(mesh2D *mesh);

void acousticsRun2D(mesh2D *mesh);

void acousticsOccaRun2D(mesh2D *mesh);

void acousticsSplitPmlOccaRun2D(mesh2D *mesh);

void acousticsVolume2D(mesh2D *mesh);
void acousticsSurface2D(mesh2D *mesh, dfloat t);
void acousticsUpdate2D(mesh2D *mesh, dfloat rka, dfloat rkb);

void acousticsPmlUpdate2D(mesh2D *mesh, dfloat rka, dfloat rkb);

void acousticsError2D(mesh2D *mesh, dfloat time);

void acousticsComputeVorticity2D(mesh2D *mesh, dfloat *q, iint outfld, iint Nfields);

void acousticsCavitySolution2D(dfloat x, dfloat y, dfloat t, dfloat *u, dfloat *v, dfloat *p);

void acousticsGaussianPulse2D(dfloat x, dfloat y, dfloat t, dfloat *u, dfloat *v, dfloat *p);

