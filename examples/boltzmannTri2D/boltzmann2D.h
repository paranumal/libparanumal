
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "mesh2D.h"


void boltzmannError2D(mesh2D *mesh, dfloat time);

void boltzmannComputeVorticity2D(mesh2D *mesh, dfloat *q, iint outfld, iint Nfields);

//dfloat boltzmannRampFunction2D(dfloat t);
void boltzmannRampFunction2D(dfloat t, dfloat *ramp, dfloat *drampdt);

// execute one Boltzmann time step using LSERK4
void boltzmannSplitPmlLserkStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer);

// execute one Boltzmann time step using LSERK4
void boltzmannSplitPmlLsimexStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer);

// execute one Boltzmann time step using LSERK4
void boltzmannSplitPmlMrabStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer);

// execute one Boltzmann time step using LSERK4
void boltzmannSplitPmlSaabStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer);

// execute one Boltzmann time step using LSERK4
void boltzmannSplitPmlSarkStep2D(mesh2D *mesh, iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer);


// output run statistics for Boltzmann simulation
void boltzmannReport2D(mesh2D *mesh, iint tstep);

//Make Perodic Connection

void boltzmannPeriodic2D(mesh2D *mesh, dfloat xper, dfloat yper);
