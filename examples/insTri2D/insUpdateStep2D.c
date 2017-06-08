#include "ins2D.h"

// complete a time step using LSERK4
void insUpdateStep2D(ins_t *ins, iint tstep, iint haloBytes,
				       dfloat * sendBuffer, dfloat * recvBuffer, 
				        char   * options){

mesh2D *mesh = ins->mesh; 

dfloat t = tstep*ins->dt + ins->dt;

//Exctract Halo On Device

if(mesh->totalHaloPairs>0){
   
     ins->totalHaloExtractKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 mesh->o_haloElementList,
                                 ins->o_PI,
                                 ins->o_pHaloBuffer);

     // copy extracted halo to HOST 
     ins->o_pHaloBuffer.copyTo(sendBuffer);            
    
     // start halo exchange
     meshHaloExchangeStart(mesh,
                           mesh->Np*sizeof(dfloat), 
                           sendBuffer,
                           recvBuffer);
    }



    // Compute Volume Contribution
   ins->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_DrT,
                            mesh->o_DsT,
                            ins->o_PI,  
                            ins->o_rhsU,
                            ins->o_rhsV);



     // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){
  // wait for halo data to arrive
    meshHaloExchangeFinish(mesh);

    ins->o_pHaloBuffer.copyFrom(recvBuffer); 

    ins->totalHaloScatterKernel(mesh->Nelements,
                                    mesh->totalHaloPairs,
                                    mesh->o_haloElementList,
                                    ins->o_PI,
                                    ins->o_pHaloBuffer);
  }




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
                              ins->PIID,
                              ins->o_PI,
                              ins->o_rhsU,
                              ins->o_rhsV);



   //computes div u^(n+1) surface term
  ins->updateUpdateKernel(mesh->Nelements,
                              ins->dt,  
                              ins->g0,
                              ins->o_U,
                              ins->o_V,
                              ins->o_P,
                              ins->o_PI,
                              ins->o_rhsU,
                              ins->o_rhsV);


   
}
