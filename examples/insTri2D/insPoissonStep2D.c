#include "ins2D.h"

// complete a time step using LSERK4
void insPoissonStep2D(ins_t *ins, iint tstep, iint haloBytes,
				       dfloat * sendBuffer, dfloat * recvBuffer, 
				        char   * options){

mesh2D *mesh = ins->mesh; 

solver_t *solver = ins->pSolver;

dfloat t = tstep*ins->dt + ins->dt;


//Exctract Halo On Device

	// if(mesh->totalHaloPairs>0){
	 
 //    ins->poissonHaloExtractKernel(mesh->Nelements,
 //                           mesh->totalHaloPairs,
 //                           mesh->o_haloElementList,
 //                           ins->o_Ux,
 //                           ins->o_Uy,
 //                           ins->o_velHaloBuffer);

 //    // copy extracted halo to HOST 
 //    ins->o_velHaloBuffer.copyTo(sendBuffer);            
 //    // start halo exchange
 //    meshHaloExchangeStart(mesh,
 //                          mesh->Np*ins->NVfields*sizeof(dfloat), 
 //                          sendBuffer,
 //                          recvBuffer);
 //  	}



   // computes div u^(n+1) volume term
   ins->divergenceVolumeKernel(mesh->Nelements,
                                 mesh->o_vgeo,
                                 mesh->o_DrT,
                                 mesh->o_DsT,
                                 ins->o_U, 
                                 ins->o_V, 
                                 ins->o_rhsP);


  //   // COMPLETE HALO EXCHANGE
  // if(mesh->totalHaloPairs>0){
  // // wait for halo data to arrive
  //   meshHaloExchangeFinish(mesh);

  //   mesh->o_haloBuffer.copyFrom(recvBuffer); 

  //   ins->poissonHaloScatterKernel(mesh->Nelements,
  //                                 mesh->totalHaloPairs,
  //                                 mesh->o_haloElementList,
  //                                 ins->o_Ux,
  //                                 ins->o_Uy,
  //                                 ins->o_velHaloBuffer);
  // }


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




  // ins->poissonRhsIpdgBCKernel(mesh->Nelements,
  //                               mesh->o_sgeo,
  //                               mesh->o_vgeo,
  //                               mesh->o_DrT,
  //                               mesh->o_DsT,
  //                               mesh->o_FMMT,
  //                               mesh->o_vmapM,
  //                               mesh->o_vmapP,
  //                               mesh->o_EToB,
  //                               t,
  //                               ins->dt,
  //                               ins->tau,
  //                               mesh->o_x,
  //                               mesh->o_y,
  //                               ins->o_Pr,
  //                               ins->o_rhsPr
  //                               );



printf("Solving for P \n");
ellipticSolveTri2D(solver, 0.0, ins->o_rhsP, ins->o_PI,  ins->pSolverOptions);

   
}
