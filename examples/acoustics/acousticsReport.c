#include "acoustics.h"

void acousticsReport(acoustics_t *acoustics, dfloat time, setupAide &newOptions){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh3D *mesh = acoustics->mesh;

  // copy data back to host
  acoustics->o_q.copyTo(mesh->q);

  // do error stuff on host
  acousticsError(mesh, time);

  // output field files
  char fname[BUFSIZ];

  sprintf(fname, "foo_%04d_%04d.vtu",rank, acoustics->frame++);

  acousticsPlotVTU(acoustics, fname);

}
