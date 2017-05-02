
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "mesh2D.h"




void boltzmannSetup2D(mesh2D *mesh, char *options);
void boltzmannRun2D(mesh2D *mesh, char *options);

void boltzmannError2D(mesh2D *mesh, dfloat time, char *opt);

void boltzmannReport2D(mesh2D *mesh, iint tstep, char *opt);
void boltzmannPeriodic2D(mesh2D *mesh, dfloat xper, dfloat yper);
void boltzmannCouetteError2D(mesh2D *mesh, dfloat time);



void boltzmannComputeVorticity2D(mesh2D *mesh, dfloat *q, iint outfld, iint Nfields);

//dfloat boltzmannRampFunction2D(dfloat t);
void boltzmannRampFunction2D(dfloat t, dfloat *ramp, dfloat *drampdt);



// Time Discretizations one step
void boltzmannLserkStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char *opt);
// Execute one Boltzmann time step using LSERK4
void boltzmannLsimexStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char *opt);

// Execute one Boltzmann time step using LSERK4
void boltzmannSark3Step2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char *opt);

// Execute one Boltzmann time step using LSERK4
void boltzmannSaab3Step2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char *opt);


