#include "ins2D.h"

// complete a time step using LSERK4
void insHelmholtzStep2D(ins_t *ins, int tstep,  int haloBytes,
				               dfloat * sendBuffer, dfloat * recvBuffer, 
                       char   * options){

	mesh2D *mesh = ins->mesh; 

  solver_t *solver = ins->vSolver; 

	dfloat t = tstep*ins->dt;


  // Exctract Halo On Device

	// if(mesh->totalHaloPairs>0){
	 
 //    ins->helmholtzHaloExtractKernel(mesh->Nelements,
 //                                    mesh->totalHaloPairs,
 //                                    mesh->o_haloElementList,
 //                                    ins->o_Ux,
 //                                    ins->o_Uy,
 //                                    ins->o_Pr,
 //                                    ins->o_totHaloBuffer);

 //    // copy extracted halo to HOST 
 //    ins->o_totHaloBuffer.copyTo(sendBuffer);            
 //    // start halo exchange
 //    meshHaloExchangeStart(mesh,
 //                          mesh->Np*(ins->NTfields)*sizeof(dfloat), // pressure also 
 //                          sendBuffer,
 //                          recvBuffer);
 //  	}

  	// Compute Volume Contribution
   ins->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_DrT,
                            mesh->o_DsT,
                            ins->o_P,  
                            ins->o_rhsU,
                            ins->o_rhsV);


  //   // COMPLETE HALO EXCHANGE
  // if(mesh->totalHaloPairs>0){
  // // wait for halo data to arrive
  //   meshHaloExchangeFinish(mesh);

  //   mesh->o_haloBuffer.copyFrom(recvBuffer); 

  //   ins->helmholtzHaloScatterKernel(mesh->Nelements,
  //                                   mesh->totalHaloPairs,
  //                                   mesh->o_haloElementList,
  //                                   ins->o_Ux,
  //                                   ins->o_Uy,
  //                                   ins->o_Pr,
  //                                   ins->o_totHaloBuffer);
  // }

 // Compute Surface Conribution
  ins->gradientSurfaceKernel(mesh->Nelements,
                              mesh->o_sgeo,
                              mesh->o_LIFTT,
                              mesh->o_vmapM,
                              mesh->o_vmapP,
                              mesh->o_EToB,
                              mesh->o_x,
                              mesh->o_y,
                              t,
                              ins->dt, 
                              ins->PID,
                              ins->o_P,
                              ins->o_rhsU,
                              ins->o_rhsV);










 // compute all forcing i.e. f^(n+1) - grad(Pr)
  ins->helmholtzRhsForcingKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_MM,
                              ins->dt,	
                              ins->a0,
                              ins->a1,
                              ins->a2,
                              ins->b0,
                              ins->b1,
                              ins->b2,
                              ins->o_U,
                              ins->o_V,
                              ins->o_UO,
                              ins->o_NU,
                              ins->o_NV,
                              ins->o_NO,
                              ins->o_rhsU,
                              ins->o_rhsV);



  
  ins->helmholtzRhsIpdgBCKernel(mesh->Nelements,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                ins->tau,
                                t,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_vgeo,
                                mesh->o_sgeo,
                                mesh->o_EToB,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                ins->o_rhsU,
                                ins->o_rhsV);

  


printf("Solving for Ux \n");
 ellipticSolveTri2D( solver, ins->lambda, ins->o_rhsU, ins->o_U, ins->vSolverOptions);

printf("Solving for Uy \n");
 ellipticSolveTri2D(solver, ins->lambda, ins->o_rhsV, ins->o_V, ins->vSolverOptions);

}
