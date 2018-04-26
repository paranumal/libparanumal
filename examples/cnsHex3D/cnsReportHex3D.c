#include "cnsHex3D.h"

void cnsReportHex3D(cns_t *cns, dfloat time, setupAide &newOptions){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh3D *mesh = cns->mesh;

  cns->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_D,
                       cns->o_q,
                       cns->o_Vort);

  // copy data back to host
  cns->o_q.copyTo(mesh->q);
  cns->o_Vort.copyTo(cns->Vort);

  // do error stuff on host
  cnsError3D(mesh, time);

  // output field files
  char fname[BUFSIZ];

  sprintf(fname, "foo_%04d_%04d.vtu",rank, cns->frame++);
  
  cnsPlotVTUHex3D(cns, fname);

}
