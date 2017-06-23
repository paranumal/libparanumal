#include "ins2D.h"

// complete a time step using LSERK4
void insPoissonStep2D(ins_t *ins, iint tstep, iint haloBytes,
				       dfloat * sendBuffer, dfloat * recvBuffer,
				        char   * options){

  mesh2D *mesh = ins->mesh;
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
                             offset,
                             ins->o_U,
                             ins->o_V,
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
                              offset,
                              ins->o_U,
                              ins->o_V,
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

  #if 0
  ins->o_rhsP.copyTo(ins->rhsP);
  dfloat maxp = 0; 
  for(iint e=0;e<mesh->Nelements;++e){
    for(iint n=0;n<mesh->Np;++n){
      const iint id = e*mesh->Np + n;
      maxp = mymax(maxp, fabs(ins->rhsP[id])); 
      const iint id2 = id+ index1*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);

      ins->P[id2] = ins->rhsP[id];
      ins->rhsP[id] = 1e-10;

    }
  }

  printf("Divergence of intermediate Velocity: %.5e \n",maxp);

 // // Just for case give exact zero divergence
  //ins->o_rhsP.copyFrom(ins->rhsP);
    
  // int Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  // ins->o_rhsP.copyTo(ins->P,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);
  #endif

  #if 1 // No time dependent BC
  ins->poissonRhsIpdgBCKernel(mesh->Nelements,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                ins->tau,
                                t,
                                ins->dt,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_vgeo,
                                mesh->o_sgeo,
                                mesh->o_EToB,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                ins->o_rhsP);
  #endif

  printf("Solving for P \n");
  //ellipticSolveTri2D(solver, 0.0, ins->o_rhsP, ins->o_PI,  ins->pSolverOptions);  
}
