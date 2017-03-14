#include "boltzmann2D.h"

void boltzmannReport2D(mesh2D *mesh, iint tstep){

  dfloat t = (tstep+1)*mesh->dt;
  
  // copy data back to host
  mesh->o_q.copyTo(mesh->q);

  // report ramp function
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0){
    dfloat ramp, drampdt;
    boltzmannRampFunction2D(t, &ramp, &drampdt);
    printf("t: %g ramp: %g drampdt: %g\n", t, ramp, drampdt);
  }
  
  // do error stuff on host
  boltzmannError2D(mesh, t);
  
  // compute vorticity
  boltzmannComputeVorticity2D(mesh, mesh->q, 0, mesh->Nfields);
  
  // output field files
  iint fld = 1;
  char fname[BUFSIZ];
  sprintf(fname, "foo_T%04d", tstep/mesh->errorStep);
  meshPlotVTU2D(mesh, fname, fld);
}
