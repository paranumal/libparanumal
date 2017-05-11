#include "ins2D.h"

// complete a time step using LSERK4
void insAdvectionStep2D(solver_t *ins,iint tstep, iint haloBytes,
				  dfloat * sendBuffer, dfloat *recvBuffer, char * options){

	mesh2D *mesh = ins->mesh; 


   // Exctract Halo On Device

	// if(mesh->totalHaloPairs>0){
 //      // EXCTRACT HALO  on DEVICE
 //      iint Nentries = mesh->Np*mesh->Nfields;
      
 //      mesh->haloExtractKernel(mesh->totalHaloPairs,
	// 		      Nentries,
	// 		      mesh->o_haloElementList,
	// 		      mesh->o_q,
	// 		      mesh->o_haloBuffer);
      
 //      // copy extracted halo to HOST 
 //      mesh->o_haloBuffer.copyTo(sendBuffer);      
      
 //      // start halo exchange
 //      meshHaloExchangeStart(mesh,
	// 		    mesh->Np*mesh->Nfields*sizeof(dfloat),
	// 		    sendBuffer,
	// 		    recvBuffer);
 //  	}

   ins->advectionVolumeKernel(
   	         mesh->Nelements,
   	         mesh->o_vgeo,
			 mesh->o_DrT,
			 mesh->o_DsT,
			 ins->o_U,
			 ins->o_rhsU
			 );


    // // COMPLETE HALO EXCHANGE
    // if(mesh->totalHaloPairs>0){
    //   // wait for halo data to arrive
    //   meshHaloExchangeFinish(mesh);
      
    //   // copy halo data to DEVICE
    //   size_t offset = mesh->Np*mesh->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
    //   mesh->o_q.copyFrom(recvBuffer, haloBytes, offset);
    // }

  

   
}
