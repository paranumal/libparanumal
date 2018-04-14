#include "cnsTri2D.h"

void cnsReportTri2D(cns_t *cns, dfloat time, char* options){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh2D *mesh = cns->mesh;

  cns->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_DrT,
                       mesh->o_DsT,
                       cns->o_q,
                       cns->o_Vort);

  // copy data back to host
  cns->o_q.copyTo(mesh->q);
  cns->o_Vort.copyTo(cns->Vort);

  // do error stuff on host
  cnsError2D(mesh, time);


  // output field files
  char fname[BUFSIZ];
  // sprintf(fname, "/u0/outputs/cns2D/foo_%04d_%04d.vtu",rank, tstep/cns->errorStep);
  sprintf(fname, "foo_%04d_%04d.vtu",rank, cns->frame++);
  //sprintf(fname, "/scratch/foo_%04d_%04d.vtu",rank, tstep/cns->errorStep);
  //  cnsPlotVTUTri2D(cns, fname);

}
