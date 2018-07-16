#include "acoustics.h"

void acousticsReport(acoustics_t *acoustics, dfloat time, setupAide &newOptions){

  mesh3D *mesh = acoustics->mesh;

  // copy data back to host
  acoustics->o_q.copyTo(mesh->q);

  // do error stuff on host
  acousticsError(mesh, time);

  // output field files
  char fname[BUFSIZ];

  sprintf(fname, "foo_%04d_%04d.vtu", mesh->rank, acoustics->frame++);

  acousticsPlotVTU(acoustics, fname);
  
}
