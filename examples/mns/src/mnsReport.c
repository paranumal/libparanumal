#include "mns.h"

void mnsReport(mns_t *mns, dfloat time, int tstep){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = mns->mesh;

  // mns->vorticityKernel(mesh->Nelements,
  //                      mesh->o_vgeo,
  //                      mesh->o_Dmatrices,
  //                      mns->fieldOffset,
  //                      mns->o_U,
  //                      mns->o_Vort);

  // mns->divergenceVolumeKernel(mesh->Nelements,
  //                            mesh->o_vgeo,
  //                            mesh->o_Dmatrices,
  //                            mns->fieldOffset,
  //                            mns->o_U,
  //                            mns->o_Div);

  // gather-scatter
  // ellipticParallelGatherScatter(mesh, mesh->ogs, ins->o_Vort, dfloatString, "add");  
  // ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_Vort, ins->o_Vort);
  
  // copy data back to host
  mns->o_U.copyTo(mns->U);
  mns->o_P.copyTo(mns->P);
  mns->o_Phi.copyTo(mns->Phi);
  // mns->o_Vort.copyTo(mns->Vort);
  // mns->o_Div.copyTo(mns->Div);

  // do error stuff on host
  mnsError(mns, time);

  if(mns->options.compareArgs("OUTPUT TYPE","VTU")){ 
    // output field files
    char fname[BUFSIZ];
    string outName;
    mns->options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(),rank, tstep/mns->outputStep);

    mnsPlotVTU(mns, fname);
  } 
}

