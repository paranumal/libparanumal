#include "ins2D.h"

// complete a time step using LSERK4
void insAdvectionStep2D(ins_t *ins, iint tstep,  iint haloBytes,
			dfloat * sendBuffer, dfloat * recvBuffer, 
			char   * options){

  mesh2D *mesh = ins->mesh; 

  dfloat t = tstep*ins->dt;


  //Exctract Halo On Device

  if(mesh->totalHaloPairs>0){
	 
     ins->totalHaloExtractKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 mesh->o_haloElementList,
                                 ins->o_U,
                                 ins->o_V,
                                 ins->o_P,
                                 ins->o_tHaloBuffer);

     // copy extracted halo to HOST 
     ins->o_tHaloBuffer.copyTo(sendBuffer);            
    
     // start halo exchange
     meshHaloExchangeStart(mesh,
                           mesh->Np*(ins->NTfields)*sizeof(dfloat), 
                           sendBuffer,
                           recvBuffer);
   	}

  // Compute Volume Contribution
  if(strstr(options, "CUBATURE")){
    ins->advectionCubatureVolumeKernel(mesh->Nelements,
				       mesh->o_vgeo,
				       mesh->o_cubDrWT,
				       mesh->o_cubDsWT,
				       mesh->o_cubInterpT,
				       ins->o_U,
				       ins->o_V,
				       ins->o_NU,
				       ins->o_NV);
  }
  else{
    ins->advectionVolumeKernel(mesh->Nelements,
			       mesh->o_vgeo,
			       mesh->o_DrT,
			       mesh->o_DsT,
			       ins->o_U,
			       ins->o_V,
			       ins->o_NU,
			       ins->o_NV);
  }

   

    // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){
  // wait for halo data to arrive
    meshHaloExchangeFinish(mesh);

    ins->o_tHaloBuffer.copyFrom(recvBuffer); 

    ins->totalHaloScatterKernel(mesh->Nelements,
                                    mesh->totalHaloPairs,
                                    mesh->o_haloElementList,
                                    ins->o_U,
                                    ins->o_V,
                                    ins->o_P,
                                    ins->o_tHaloBuffer);
  }

  // Compute Volume Contribution
  if(strstr(options, "CUBATURE")){
    // Compute Surface Conribution
    ins->advectionCubatureSurfaceKernel(mesh->Nelements,
					mesh->o_sgeo,
					mesh->o_intInterpT,
					mesh->o_intLIFTT,
					mesh->o_vmapM,
					mesh->o_vmapP,
					mesh->o_EToB,
					t,
					mesh->o_intx,
					mesh->o_inty,
					ins->o_U,
					ins->o_V,
					ins->o_NU,
					ins->o_NV);
  }
  else{
    // Compute Surface Conribution
    ins->advectionSurfaceKernel(mesh->Nelements,
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
				ins->o_NU,
				ins->o_NV);
  }





 
}
