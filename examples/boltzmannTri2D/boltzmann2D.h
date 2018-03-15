
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "mesh2D.h"
#include <complex.h>  



void boltzmannSetup2D(mesh2D *mesh, char *options);
//void boltzmannSetupTest2D(mesh2D *mesh, char *options);

void boltzmannMRABSetup2D(mesh2D *mesh, char *options);


void boltzmannMRABPmlSetup2D(mesh2D *mesh, char *options);
void boltzmannPmlSetup2D(mesh2D *mesh, char *options);


void boltzmannRun2D(mesh2D *mesh, char *options);
// void boltzmannRunTest2D(mesh2D *mesh, char *options);

// void boltzmannMRRun2D(mesh2D *mesh, char *options);

void boltzmannError2D(mesh2D *mesh, dfloat time, char *opt);
void boltzmannForces2D(mesh2D *mesh, dfloat time, char *opt);

void boltzmannReport2D(mesh2D *mesh, int tstep, char *opt);
void boltzmannPeriodic2D(mesh2D *mesh, dfloat xper, dfloat yper);
void boltzmannCouetteError2D(mesh2D *mesh, dfloat time);

void boltzmannPlotVTUField2D(mesh2D *mesh, char *fname);
void boltzmannPlotVTU2D(mesh2D *mesh, char * FileName);
void boltzmannPlotTEC2D(mesh2D *mesh, char * FileName, dfloat solutionTime);
void boltzmannPlotFieldTEC2D(mesh2D *mesh, char * FileName, dfloat solutionTime, int fld);
void boltzmannComputeVorticity2D(mesh2D *mesh, dfloat *q, int outfld, int Nfields);

//dfloat boltzmannRampFunction2D(dfloat t);
void boltzmannRampFunction2D(dfloat t, dfloat *ramp, dfloat *drampdt);

//
//void boltzmannPulse2D(mesh2D *mesh, int id, int cnt);


void boltzmannSAABStep2D(mesh2D *mesh, int tstep, int haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char *opt);


// / Time Discretizations one step LSERK
void boltzmannSRABStep2D(mesh2D *mesh, int tstep, int haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char *opt);

// Time Discretizations one step LSERK
void boltzmannLSERKStep2D(mesh2D *mesh, int tstep, int haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char *opt);
// Execute one Boltzmann time step using LSIMEX
void boltzmannLSIMEXStep2D(mesh2D *mesh, int tstep, int haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char *opt);

//Execute one Boltzmann time step using SARK
void boltzmannSARKStep2D(mesh2D *mesh, int tstep, int haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char *opt);

// // Execute one Boltzmann time step using SAAB
// void boltzmannSaab3Step2D(mesh2D *mesh, int tstep, int haloBytes,
// 				  dfloat * sendBuffer, dfloat *recvBuffer, char *opt);

// Execute one Boltzmann time step using MRAB
void boltzmannMRABStep2D(mesh2D *mesh, int tstep, int haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char *opt);

// Execute one Boltzmann time step using MRSAAB
void boltzmannMRSAABStep2D(mesh2D *mesh, int tstep, int haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char *opt);



