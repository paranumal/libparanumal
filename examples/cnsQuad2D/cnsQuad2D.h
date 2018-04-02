#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"

#include "mesh2D.h"

void cnsRunQuad2D(mesh2D *mesh);

void cnsSetupQuad2D(mesh2D *mesh);

void cnsVolumeQuad2D(mesh2D *mesh);

void cnsSurfaceQuad2D(mesh2D *mesh);

void cnsUpdate2D(mesh2D *mesh, dfloat rka, dfloat rkb);

void cnsError2D(mesh2D *mesh, dfloat time);

void cnsCavitySolution2D(dfloat x, dfloat y, dfloat t,
			 dfloat *u, dfloat *v, dfloat *p);


void cnsGaussianPulse2D(dfloat x, dfloat y, dfloat t,
			 dfloat *u, dfloat *v, dfloat *p);
