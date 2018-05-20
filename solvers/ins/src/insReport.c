#include "ins.h"

void insReport(ins_t *ins, dfloat time, int tstep){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = ins->mesh;

  ins->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_Dmatrices,
                       ins->fieldOffset,
                       ins->o_U,
                       ins->o_Vort);

  ins->divergenceVolumeKernel(mesh->Nelements,
                             mesh->o_vgeo,
                             mesh->o_Dmatrices,
                             ins->fieldOffset,
                             ins->o_U,
                             ins->o_Div);

  // gather-scatter
  ellipticParallelGatherScatter(mesh, mesh->ogs, ins->o_Vort, dfloatString, "add");  
  ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_Vort, ins->o_Vort);
  
  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_P.copyTo(ins->P);

  ins->o_Vort.copyTo(ins->Vort);
  ins->o_Div.copyTo(ins->Div);

  // do error stuff on host
  insError(ins, time);

  if(ins->options.compareArgs("OUTPUT TYPE","VTU")){ 
    // output field files
    char fname[BUFSIZ];
    string outName;
    ins->options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d_%04d.vtu",(char*)outName.c_str(),rank, tstep/ins->outputStep);

    insPlotVTU(ins, fname);
  } 
}

