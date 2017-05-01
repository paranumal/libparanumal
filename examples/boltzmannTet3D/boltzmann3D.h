#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "mesh3D.h"



void boltzmannSetup3D(mesh3D *mesh, char *opt);

void boltzmannRun3D(mesh3D *mesh, char *opt);

void boltzmannReport3D(mesh3D *mesh, iint tstep, char *options);

void boltzmannError3D(mesh3D *mesh, dfloat time,char *options);

void boltzmannRampFunction3D(dfloat t, dfloat *ramp, dfloat *drampdt);

void boltzmannLserkStep3D(mesh3D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char *opt);