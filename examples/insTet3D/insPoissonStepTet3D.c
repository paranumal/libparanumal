#include "insTet3D.h"

// complete a time step using LSERK4
void insPoissonStepTet3D(ins_t *ins, int tstep, const char* options){

  mesh3D *mesh = ins->mesh;
  solver_t *solver = ins->pSolver;
  dfloat t = tstep*ins->dt + ins->dt;

  //hard coded for 3 stages.
  //The result of the helmholtz solve is stored in the next index
  int index1 = (ins->index+1)%3;
  dlong offset  = mesh->Nelements+mesh->totalHaloPairs;
  dlong ioffset = index1*offset;

  if(mesh->totalHaloPairs>0){
    ins->velocityHaloExtractKernel(mesh->Nelements,
                               mesh->totalHaloPairs,
                               mesh->o_haloElementList,
                               offset,
                               ins->o_U,
                               ins->o_V,
                               ins->o_W,
                               ins->o_vHaloBuffer);

    // copy extracted halo to HOST 
    ins->o_vHaloBuffer.copyTo(ins->vSendBuffer);           
  
    // start halo exchange
    meshHaloExchangeStart(mesh,
                         mesh->Np*(ins->NVfields)*sizeof(dfloat),
                         ins->vSendBuffer,
                         ins->vRecvBuffer);
  }

  // computes div u^(n+1) volume term
  ins->divergenceVolumeKernel(mesh->Nelements,
                             mesh->o_vgeo,
                             mesh->o_DrT,
                             mesh->o_DsT,
                             mesh->o_DtT,
                             offset,
                             ins->o_U,
                             ins->o_V,
                             ins->o_W,
                             ins->o_rhsP);

  if(mesh->totalHaloPairs>0){
    meshHaloExchangeFinish(mesh);

    ins->o_vHaloBuffer.copyFrom(ins->vRecvBuffer); 

    ins->velocityHaloScatterKernel(mesh->Nelements,
                                  mesh->totalHaloPairs,
                                  mesh->o_haloElementList,
                                  offset,
                                  ins->o_U,
                                  ins->o_V,
                                  ins->o_W,
                                  ins->o_vHaloBuffer);
  }


  //computes div u^(n+1) surface term
  ins->divergenceSurfaceKernel(mesh->Nelements,
                              mesh->o_sgeo,
                              mesh->o_LIFTT,
                              mesh->o_vmapM,
                              mesh->o_vmapP,
                              mesh->o_EToB,
                              t,
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              offset,
                              ins->o_U,
                              ins->o_V,
                              ins->o_W,
                              ins->o_rhsP);

  // compute all forcing i.e. f^(n+1) - grad(Pr)
  ins->poissonRhsForcingKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_MM,
                              ins->dt,  
                              ins->g0,
                              ins->o_rhsP);

#if 0
  //add penalty from jumps in previous pressure
  ins->poissonPenaltyKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_vgeo,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_DtT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                ins->tau,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                t,
                                ins->dt,
                                ins->c0,
                                ins->c1,
                                ins->c2,
                                ins->index,
                                (mesh->Nelements+mesh->totalHaloPairs),
                                ins->o_P,
                                ins->o_rhsP);
  #endif

  #if 1// if time dependent BC or Pressure Solve not Increment
  const int pressure_solve = 0; // ALGEBRAIC SPLITTING 
  if (strstr(ins->pSolverOptions,"CONTINUOUS")) {
    ins->poissonRhsBCKernel(mesh->Nelements,
                            pressure_solve,
                            mesh->o_ggeo,
                            mesh->o_sgeo,
                            mesh->o_SrrT,
                            mesh->o_SrsT,
                            mesh->o_SrtT,
                            mesh->o_SsrT,
                            mesh->o_SssT,
                            mesh->o_SstT,
                            mesh->o_StrT,
                            mesh->o_StsT,
                            mesh->o_SttT,
                            mesh->o_MM,
                            mesh->o_vmapM,
                            mesh->o_sMT,
                            t,
                            ins->dt,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            ins->o_PmapB,
                            ins->o_rhsP);

    // gather-scatter
    ellipticParallelGatherScatterTet3D(mesh, mesh->ogs, ins->o_rhsP, dfloatString, "add");  
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_rhsP);

  } else if (strstr(ins->pSolverOptions,"IPDG")) {
    ins->poissonRhsIpdgBCKernel(mesh->Nelements,
                                pressure_solve,
                                mesh->o_vmapM,
                                ins->tau,
                                t,
                                ins->dt,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                mesh->o_vgeo,
                                mesh->o_sgeo,
                                mesh->o_EToB,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_DtT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                ins->o_rhsP);
  }
  #endif

  ins->NiterP = ellipticSolveTet3D(solver, 0.0, ins->presTOL, ins->o_rhsP, ins->o_PI,  ins->pSolverOptions); 

  if (strstr(ins->pSolverOptions,"CONTINUOUS")) {
    ins->poissonAddBCKernel(mesh->Nelements,
                            pressure_solve,
                            t,
                            ins->dt,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            ins->o_PmapB,
                            ins->o_PI);
  }
}
