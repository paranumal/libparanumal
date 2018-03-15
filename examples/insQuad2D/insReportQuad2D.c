#include "insQuad2D.h"

void insReportQuad2D(ins_t *ins, int tstep, char *options){
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh2D *mesh = ins->mesh;

  dfloat t = (tstep)*ins->dt;
  
  dlong offset = ins->index*(mesh->Nelements+mesh->totalHaloPairs);
  ins->vorticityKernel(mesh->Nelements,
                       mesh->o_vgeo,
                       mesh->o_D,
                       offset,
                       ins->o_U,
                       ins->o_V,
                       ins->o_Vort);

  ins->divergenceVolumeKernel(mesh->Nelements,
                             mesh->o_vgeo,
                             mesh->o_D,
                             offset,
                             ins->o_U,
                             ins->o_V,
                             ins->o_Div);

  // gather-scatter
  ellipticParallelGatherScatterQuad2D(mesh, mesh->ogs, ins->o_Vort, dfloatString, "add");  
  ins->pSolver->dotMultiplyKernel(mesh->Nelements*mesh->Np, mesh->ogs->o_invDegree, ins->o_Vort, ins->o_Vort);

  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_V.copyTo(ins->V); 
  ins->o_P.copyTo(ins->P);

  ins->o_Vort.copyTo(ins->Vort);
  ins->o_Div.copyTo(ins->Div);
  
  // do error stuff on host
  insErrorQuad2D(ins, t, options);
 
  if(strstr(options, "VTU")){ 
    // output field files
    char fname[BUFSIZ];
    // sprintf(fname, "/u0/outputs/ins2D/foo_%04d_%04d.vtu",rank, tstep/ins->errorStep);
    sprintf(fname, "foo_%04d_%04d.vtu",rank, tstep/ins->errorStep);

    insPlotVTUQuad2D(ins, fname);
  } 
}

