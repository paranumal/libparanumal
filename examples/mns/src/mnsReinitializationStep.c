#include "mns.h"

// complete a time step using LSERK4
void mnsReinitializationStep(mns_t *mns, int tstep, int haloBytes, dfloat * sendBuffer, dfloat *recvBuffer){



mesh_t *mesh = mns->mesh; 


// field offset at this step // index = 0 for levelSet run only, 
dlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

// // LSERK4 stages
// for(int rk=0;rk<mesh->Nrk;++rk){

//   // intermediate stage time
//   dfloat time = tstep*mns->dt + mns->dt*mesh->rkc[rk];

//   if(mesh->totalHaloPairs>0){
//     #if ASYNC 
//       mesh->device.setStream(dataStream);
//     #endif

//     mesh->haloExtractKernel(mesh->totalHaloPairs,
//                 mesh->Np,
//                 mesh->o_haloElementList,
//                 mns->o_Phi,
//                 mesh->o_haloBuffer);

//     // copy extracted halo to HOST
//     mesh->o_haloBuffer.asyncCopyTo(sendBuffer);

//     #if ASYNC 
//       mesh->device.setStream(defaultStream);
//     #endif
//   }
    
//     occaTimerTic(mesh->device, "LevelSetVolumeKernel");   
//     mns->levelSetVolumeKernel(mesh->Nelements,
//                               mesh->o_vgeo,
//                               mesh->o_cubvgeo,
//                               mesh->o_cubDWmatrices,
//                               mesh->o_cubInterpT,
//                               mesh->o_cubProjectT,
//                               mns->fieldOffset,
//                               time, 
//                               mns->o_U,
//                               mns->o_Phi,
//                               mns->o_rhsPhi);         

//     occaTimerToc(mesh->device, "LevelSetVolumeKernel");   

//      if(mesh->totalHaloPairs>0){
    
//     #if ASYNC 
//       mesh->device.setStream(dataStream);
//     #endif

//     //make sure the async copy is finished
//     mesh->device.finish();
//     // start halo exchange
//     meshHaloExchangeStart(mesh,
//                           mesh->Np*sizeof(dfloat),
//                           sendBuffer,
//                           recvBuffer);
//     // wait for halo data to arrive
//     meshHaloExchangeFinish(mesh);
//     // copy halo data to DEVICE
//     size_t offset = mesh->Np*mesh->Nelements*sizeof(dfloat); // offset for halo data
//     mns->o_Phi.asyncCopyFrom(recvBuffer, haloBytes, offset);
//     mesh->device.finish();        

//     #if ASYNC 
//       mesh->device.setStream(defaultStream);
//     #endif
//   }


//       occaTimerTic(mesh->device,"LevelSetSurfaceKernel");
//       mns->levelSetSurfaceKernel(mesh->Nelements,
//                                 mesh->o_vgeo,
//                                 mesh->o_sgeo,
//                                 mesh->o_cubsgeo,
//                                 mesh->o_intInterpT,
//                                 mesh->o_intLIFTT,
//                                 mesh->o_cubInterpT,
//                                 mesh->o_cubProjectT,
//                                 mesh->o_vmapM,
//                                 mesh->o_vmapP,
//                                 mesh->o_EToB,
//                                 time,
//                                 mesh->o_intx,
//                                 mesh->o_inty,
//                                 mesh->o_intz,
//                                 mns->fieldOffset,
//                                 mns->o_U,
//                                 mns->o_Phi,
//                                 mns->o_rhsPhi);  
      
//       occaTimerToc(mesh->device,"LevelSetSurfaceKernel");
 
 

 
//     occaTimerTic(mesh->device,"LevelSetUpdateKernel");
//     mns->levelSetUpdateKernel(mesh->Nelements,
//                               mns->dt,
//                               mesh->rka[rk],
//                               mesh->rkb[rk],
//                               mns->o_rhsPhi,
//                               mns->o_resPhi,
//                               mns->o_Phi);
//     occaTimerToc(mesh->device,"LevelSetUpdateKernel");
//  }


}