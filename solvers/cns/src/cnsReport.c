#include "cns.h"

void cnsReport(cns_t *cns, dfloat time, setupAide &options){

  mesh3D *mesh = cns->mesh;

  cns->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_Dmatrices,
                       cns->o_q,
                       cns->o_Vort);

  // copy data back to host
  cns->o_q.copyTo(mesh->q);
  cns->o_Vort.copyTo(cns->Vort);

  // do error stuff on host
  cnsError(mesh, time);

  //  cnsForces(cns, time);

  // output field files
  char fname[BUFSIZ];
  string outName;
  options.getArgs("OUTPUT FILE NAME", outName);
  sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(), mesh->rank, cns->frame++);

  cnsPlotVTU(cns, fname);

}
