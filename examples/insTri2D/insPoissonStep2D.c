#include "ins2D.h"

// complete a time step using LSERK4
void insPoissonStep2D(ins_t *ins, iint tstep, iint haloBytes,
				       dfloat * sendBuffer, dfloat * recvBuffer, 
				        char   * options){

mesh2D *mesh = ins->mesh; 

solver_t *solver = ins->prsolver;

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
                                 ins->o_Ux, 
                                 ins->o_Uy, 
                                 ins->o_rhsPr);


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
                              ins->o_Ux,
                              ins->o_Uy,
                              ins->o_rhsPr);




  // compute all forcing i.e. f^(n+1) - grad(Pr)
  ins->poissonRhsForcingKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_MM,
                              ins->dt,  
                              ins->g0,
                              ins->o_rhsPr);




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



// dfloat *Pr2    = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*ins->Nfields,sizeof(dfloat));

// ins->o_rhsPr.copyTo(ins->Pr);
//o_PrI.copyFrom(Pr2);
ellipticSolveTri2D(solver, 0.0, ins->o_rhsPr, ins->o_PrI,  ins->prsolverOptions);

   
}
