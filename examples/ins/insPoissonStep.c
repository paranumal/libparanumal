#include "ins.h"

// complete a time step using LSERK4
void insPoissonStep(ins_t *ins, dfloat time){

  mesh2D *mesh = ins->mesh;
  solver_t *solver = ins->pSolver;
  dfloat time = tstep*ins->dt + ins->dt;

  dlong offset  = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;

  /* note: the surface kernel isn't needed with continuous pressure. Just the inflow boundary 
           contributions to the surface 
           TODO: Need a separate kernel to do the surface kernel for just boundaries */
  //if (strstr(ins->pSolverOptions,"IPDG")) {
    if(mesh->totalHaloPairs>0){
      ins->velocityHaloExtractKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 mesh->o_haloElementList,
                                 offset,
                                 ins->o_U,
                                 ins->o_vHaloBuffer);

      // copy extracted halo to HOST 
      ins->o_vHaloBuffer.copyTo(ins->vSendBuffer);           
    
      // start halo exchange
      meshHaloExchangeStart(mesh,
                           mesh->Np*(ins->NVfields)*sizeof(dfloat),
                           ins->vSendBuffer,
                           ins->vRecvBuffer);
    }
  //}
  
  occaTimerTic(mesh->device,"DivergenceVolume");
  // computes div u^(n+1) volume term
  ins->divergenceVolumeKernel(mesh->Nelements,
                             mesh->o_vgeo,
                             mesh->o_Dmatrices,
                             offset,
                             ins->o_U,
                             ins->o_rhsP);
   occaTimerToc(mesh->device,"DivergenceVolume");

  //if (strstr(ins->pSolverOptions,"IPDG")) {
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);

      ins->o_vHaloBuffer.copyFrom(ins->vRecvBuffer); 

      ins->velocityHaloScatterKernel(mesh->Nelements,
                                    mesh->totalHaloPairs,
                                    offset,
                                    ins->o_U,
                                    ins->o_vHaloBuffer);
    }

    occaTimerTic(mesh->device,"DivergenceSurface");
    //computes div u^(n+1) surface term
    ins->divergenceSurfaceKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_LIFTT,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                time,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                offset,
                                ins->o_U,
                                ins->o_rhsP);
    occaTimerToc(mesh->device,"DivergenceSurface");
  //}

  
  occaTimerTic(mesh->device,"PoissonRhsForcing");
  // compute all forcing i.e. f^(n+1) - grad(Pr)
  ins->poissonRhsForcingKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_MM,
                              ins->dt,  
                              ins->g0,
                              ins->o_rhsP);
  occaTimerToc(mesh->device,"PoissonRhsForcing");

#if 0
  //add penalty from jumps in previous pressure
  ins->poissonPenaltyKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_vgeo,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                ins->tau,
                                mesh->o_x,
                                mesh->o_y,
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

  #if 1 // if time dependent BC
  //
  const int pressure_solve = 0; // ALGEBRAIC SPLITTING 
  if (strstr(ins->pSolverOptions,"CONTINUOUS")) {
    ins->poissonRhsBCKernel(mesh->Nelements,
                            pressure_solve,
                            mesh->o_ggeo,
                            mesh->o_sgeo,
                            mesh->o_Dmatrices,
                            mesh->o_Smatrices,
                            mesh->o_MM,
                            mesh->o_vmapM,
                            mesh->o_sMT,
                            time,
                            ins->dt,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            ins->o_PmapB,
                            ins->o_rhsP);
  } else if (strstr(ins->pSolverOptions,"IPDG")) {
    occaTimerTic(mesh->device,"PoissonRhsIpdg"); 
    ins->poissonRhsIpdgBCKernel(mesh->Nelements,
                                  pressure_solve,
                                  mesh->o_vmapM,
                                  ins->tau,
                                  time,
                                  ins->dt,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_Dmatrices,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  ins->o_rhsP);
    occaTimerToc(mesh->device,"PoissonRhsIpdg");
  }
  #endif

  //keep current PI as the initial guess?

  // gather-scatter
  if(strstr(ins->pSolverOptions, "CONTINUOUS")){
    ellipticParallelGatherScatter(mesh, mesh->ogs, ins->o_rhsP, dfloatString, "add");  
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_rhsP);
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_PI);
  }

  occaTimerTic(mesh->device,"Pr Solve");
  ins->NiterP = ellipticSolve(solver, 0.0, ins->presTOL, ins->o_rhsP, ins->o_PI); 
  occaTimerToc(mesh->device,"Pr Solve"); 

  if (strstr(ins->pSolverOptions,"CONTINUOUS")) {
    ins->poissonAddBCKernel(mesh->Nelements,
                            pressure_solve,
                            time,
                            ins->dt,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            ins->o_PmapB,
                            ins->o_PI);
  }
}
