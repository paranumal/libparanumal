#include "insHex3D.h"

void insReportHex3D(ins_t *ins, int tstep, char *options){
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh3D *mesh = ins->mesh;

  dfloat t = (tstep)*ins->dt;
  
  dlong offset = ins->index*(mesh->Nelements+mesh->totalHaloPairs);
  ins->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_D,
                       offset,
                       ins->o_U,
                       ins->o_V,
                       ins->o_W,
                       ins->o_Vx,
                       ins->o_Vy,
                       ins->o_Vz);

  ins->divergenceVolumeKernel(mesh->Nelements,
                             mesh->o_vgeo,
                             mesh->o_D,
                             offset,
                             ins->o_U,
                             ins->o_V,
                             ins->o_W,
                             ins->o_Div);

  // gather-scatter
  ellipticParallelGatherScatterHex3D(mesh, mesh->ogs, ins->o_Vx, dfloatString, "add");
  ellipticParallelGatherScatterHex3D(mesh, mesh->ogs, ins->o_Vy, dfloatString, "add");  
  ellipticParallelGatherScatterHex3D(mesh, mesh->ogs, ins->o_Vz, dfloatString, "add");  
  ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_Vx, ins->o_Vx);
  ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_Vy, ins->o_Vy);
  ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_Vz, ins->o_Vz);  

  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_V.copyTo(ins->V); 
  ins->o_V.copyTo(ins->W); 
  ins->o_P.copyTo(ins->P);

  ins->o_Vx.copyTo(ins->Vx);
  ins->o_Vy.copyTo(ins->Vy);
  ins->o_Vz.copyTo(ins->Vz);
  ins->o_Div.copyTo(ins->Div);
  
  // do error stuff on host
  insErrorHex3D(ins, t, options);
 
  if(strstr(options, "VTU")){ 
    // output field files
    char fname[BUFSIZ];
    // sprintf(fname, "/u0/outputs/ins3D/foo_%04d_%04d.vtu",rank, tstep/ins->errorStep);
    sprintf(fname, "foo_%04d_%04d.vtu",rank, tstep/ins->errorStep);

    insPlotVTUHex3D(ins, fname);
  } 
}

