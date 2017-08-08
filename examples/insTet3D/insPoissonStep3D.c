#include "ins3D.h"

// complete a time step using LSERK4
void insPoissonStep3D(ins_t *ins, iint tstep, iint haloBytes,
				       dfloat * sendBuffer, dfloat * recvBuffer,
				        char   * options){

  mesh3D *mesh = ins->mesh;
  solver_t *solver = ins->pSolver;
  dfloat t = tstep*ins->dt + ins->dt;

  //hard coded for 3 stages.
  //The result of the helmholtz solve is stored in the next index
  int index1 = (ins->index+1)%3;

  iint offset = index1*(mesh->Nelements+mesh->totalHaloPairs);

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
    ins->o_vHaloBuffer.copyTo(sendBuffer);           
  
    // start halo exchange
    meshHaloExchangeStart(mesh,
                         mesh->Np*(ins->NVfields)*sizeof(dfloat),
                         sendBuffer,
                         recvBuffer);
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

    ins->o_vHaloBuffer.copyFrom(recvBuffer); 

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

  #if 0// if time dependent BC or Pressure Solve not Increment
  ins->poissonRhsIpdgBCKernel(mesh->Nelements,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
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
  #endif

  //ins->o_rhsP.copyTo(ins->rhsP);

  // // 
  // dfloat maxrhsp = 0; 
  // for (iint e=0; e<mesh->Nelements; e++){
  //   for (iint n=0; n<mesh->Np; n++){
  //     maxrhsp = mymax(maxrhsp, ins->rhsP[n+e*mesh->Np]);
  //     //ins->rhsP[n+e*mesh->Np] = 0.00000001; 
  //   }
  // }
  
  // printf("maxRhsP = %.5e, dt= %.5e\n",maxrhsp, ins->dt);
  // ins->o_rhsP.copyFrom(ins->rhsP);

  printf("Solving for P \n");
  ellipticSolveTet3D(solver, 0.0, ins->o_rhsP, ins->o_PI,  ins->pSolverOptions); 

  //ins->o_PI.copyFrom(ins->rhsP); 
}
