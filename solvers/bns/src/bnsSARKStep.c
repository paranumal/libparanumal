/*

  The MIT License (MIT)

  Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

*/

#include "bns.h"

#define BNS_ANSYC 0

// complete a time step using LSERK4
void bnsSARKStep(bns_t *bns, dfloat time, int haloBytes,
		 dfloat * sendBuffer, dfloat *recvBuffer, setupAide &options){


  // bns->shiftIndex = 0; 
  mesh_t *mesh = bns->mesh; 

  dlong offset    = mesh->Nelements*mesh->Np*bns->Nfields;
  dlong pmloffset = mesh->pmlNelements*mesh->Np*bns->Nfields;

  const dlong  dzero = 0.0; 
  const int    izero = 0; 

  // LSERK4 stages
  for(int rk=0;rk<bns->NrkStages;++rk){

    // intermediate stage time
    dfloat currentTime = time + bns->rkC[rk]*bns->dt;

    occaTimerTic(mesh->device, "RKStageKernel");  
    if(mesh->nonPmlNelements){
      occaTimerTic(mesh->device, "NonPmlRKStageKernel");  
      bns->updateStageKernel(mesh->nonPmlNelements,
			     mesh->o_nonPmlElementIds,
			     offset,
			     rk,
			     bns->dt,
			     bns->o_sarkC,
			     bns->o_rkA,
			     bns->o_sarkA,
			     bns->o_q,
			     bns->o_rkrhsq,
			     bns->o_rkq);
      occaTimerToc(mesh->device, "NonPmlRKStageKernel");  
    }
  
    if(mesh->pmlNelements){
      occaTimerTic(mesh->device, "PmlRKStageKernel");  
      bns->pmlUpdateStageKernel(mesh->pmlNelements,
                                mesh->o_pmlElementIds,
                                mesh->o_pmlIds,
                                offset,
                                pmloffset,
                                rk,
                                bns->dt,
                                bns->o_sarkC,
                                bns->o_rkA,
                                bns->o_sarkA,
                                bns->o_q,
                                bns->o_pmlqx,
                                bns->o_pmlqy,
                                bns->o_pmlqz,
                                bns->o_rkrhsq,
                                bns->o_rkrhsqx,
                                bns->o_rkrhsqy,
                                bns->o_rkrhsqz,
                                bns->o_rkq,
                                bns->o_rkqx,
                                bns->o_rkqy,
                                bns->o_rkqz);
      occaTimerToc(mesh->device, "PmlRKStageKernel");
    }

    occaTimerToc(mesh->device, "RKStageKernel");  

    if(mesh->totalHaloPairs>0){
#if BNS_ASYNC 
      mesh->device.setStream(dataStream);

      int Nentries = mesh->Np*bns->Nfields;
      mesh->haloExtractKernel(mesh->totalHaloPairs,
			      Nentries,
			      mesh->o_haloElementList,
			      bns->o_rkq,
			      mesh->o_haloBuffer);

      // copy extracted halo to HOST
      mesh->o_haloBuffer.copyTo(sendBuffer,"async: true");
      mesh->device.setStream(defaultStream);

#else
      int Nentries = mesh->Np*bns->Nfields;
      mesh->haloExtractKernel(mesh->totalHaloPairs,
			      Nentries,
			      mesh->o_haloElementList,
			      bns->o_rkq,
			      mesh->o_haloBuffer);
      // copy extracted halo to HOST
      mesh->o_haloBuffer.copyTo(sendBuffer);
      // start halo exchange
      meshHaloExchangeStart(mesh, bns->Nfields*mesh->Np*sizeof(dfloat), sendBuffer, recvBuffer);
#endif
    }
    
    // dfloat ramp = 1.0, drampdt = 0.0; 
    // COMPUTE RAMP FUNCTION 
    dfloat fx, fy, fz, intfx, intfy, intfz;
    bnsBodyForce(currentTime , &fx, &fy, &fz, &intfx, &intfy, &intfz);

    occaTimerTic(mesh->device, "VolumeKernel");    
    // compute volume contribution to DG boltzmann RHS
    if(mesh->pmlNelements){ 
      occaTimerTic(mesh->device,"PmlVolumeKernel");

      if(bns->pmlcubature){
        bns->pmlVolumeKernel(mesh->pmlNelements,
                             mesh->o_pmlElementIds,
                             mesh->o_pmlIds,
                             dzero,
                             dzero,
                             izero,
                             fx,fy,fz,
                             mesh->o_vgeo,
                             mesh->o_Dmatrices,
                             bns->o_rkq,
                             bns->o_rkqx,
                             bns->o_rkqy,
                             bns->o_rkqz,
                             bns->o_rhsq,
                             bns->o_pmlrhsqx,
                             bns->o_pmlrhsqy,
                             bns->o_pmlrhsqz);
      }else{
        bns->pmlVolumeKernel(mesh->pmlNelements,
                             mesh->o_pmlElementIds,
                             mesh->o_pmlIds,
                             dzero,
                             dzero,
                             izero,
                             fx,fy,fz,
                             mesh->o_vgeo,
                             bns->o_pmlSigmaX,
                             bns->o_pmlSigmaY,
                             bns->o_pmlSigmaZ,
                             mesh->o_Dmatrices,
                             bns->o_rkq,
                             bns->o_rkqx,
                             bns->o_rkqy,
                             bns->o_rkqz,
                             bns->o_rhsq,
                             bns->o_pmlrhsqx,
                             bns->o_pmlrhsqy,
                             bns->o_pmlrhsqz);
      }
      
      occaTimerToc(mesh->device,"PmlVolumeKernel");

    }

    // compute volume contribution to DG boltzmann RHS added d/dt (ramp(qbar)) to RHS
    if(mesh->nonPmlNelements){
      occaTimerTic(mesh->device,"NonPmlVolumeKernel");
      bns->volumeKernel(mesh->nonPmlNelements,
			mesh->o_nonPmlElementIds,
			dzero,
			izero,
			fx,fy, fz,
			mesh->o_vgeo,
			mesh->o_x,
			mesh->o_y,
			mesh->o_z,
			mesh->o_Dmatrices,
			bns->o_rkq,
			bns->o_rhsq);
      occaTimerToc(mesh->device,"NonPmlVolumeKernel");
    }
    occaTimerToc(mesh->device, "VolumeKernel");    
    

#if 1

    occaTimerTic(mesh->device, "RelaxationKernel");
    if(mesh->pmlNelements){
      occaTimerTic(mesh->device, "PmlRelaxationKernel");

      if(bns->pmlcubature){
	bns->pmlRelaxationKernel(mesh->pmlNelements,
				 mesh->o_pmlElementIds,
				 mesh->o_pmlIds,
				 mesh->o_vgeo,
				 mesh->o_cubvgeo,
				 dzero,
				 dzero,
				 izero,
				 mesh->o_cubInterpT,
				 mesh->o_cubProjectT,
				 bns->o_pmlSigmaX,
				 bns->o_pmlSigmaY,
				 bns->o_pmlSigmaZ,
				 bns->o_rkq,
				 bns->o_rkqx,
				 bns->o_rkqy,
				 bns->o_rkqz,
				 bns->o_rhsq,
				 bns->o_pmlrhsqx,
				 bns->o_pmlrhsqy,
				 bns->o_pmlrhsqz);  
      }else{
	bns->pmlRelaxationKernel(mesh->pmlNelements,
				 mesh->o_pmlElementIds,
				 mesh->o_pmlIds,
				 mesh->o_vgeo,
				 mesh->o_cubvgeo,
				 dzero,
				 dzero,
				 izero,
				 mesh->o_cubInterpT,
				 mesh->o_cubProjectT,
				 bns->o_rkq,
				 bns->o_rhsq);

      }
     
      occaTimerToc(mesh->device, "PmlRelaxationKernel");
    }

    // compute relaxation terms using cubature
    if(mesh->nonPmlNelements){
      occaTimerTic(mesh->device, "NonPmlRelaxationKernel");
      bns->relaxationKernel(mesh->nonPmlNelements,
			    mesh->o_nonPmlElementIds,
			    mesh->o_vgeo,
			    mesh->o_cubvgeo,
			    dzero, // 0
			    izero, //0
			    mesh->o_cubInterpT,
			    mesh->o_cubProjectT,
			    bns->o_rkq,
			    bns->o_rhsq);  
      occaTimerToc(mesh->device, "NonPmlRelaxationKernel");
    }
    // VOLUME KERNELS
    occaTimerToc(mesh->device, "RelaxationKernel");
#endif
    
    if(mesh->totalHaloPairs>0){
    
#if BNS_ASYNC 
      mesh->device.setStream(dataStream);

      //make sure the async copy is finished
      mesh->device.finish();
      // start halo exchange
      meshHaloExchangeStart(mesh,
			    bns->Nfields*mesh->Np*sizeof(dfloat),
			    sendBuffer,
			    recvBuffer);
      // wait for halo data to arrive
      meshHaloExchangeFinish(mesh);
      // copy halo data to DEVICE
      size_t offset = mesh->Np*bns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      bns->o_rkq.copyFrom(recvBuffer, haloBytes, offset,"async: true");
      mesh->device.finish();

      mesh->device.setStream(defaultStream);
#else
      meshHaloExchangeFinish(mesh);
      // copy halo data to DEVICE
      size_t offset = mesh->Np*bns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      bns->o_rkq.copyFrom(recvBuffer, haloBytes, offset);
#endif

    }



    // SURFACE KERNELS
    occaTimerTic(mesh->device,"SurfaceKernel");

    if(mesh->pmlNelements){
      occaTimerTic(mesh->device,"PmlSurfaceKernel");
      bns->pmlSurfaceKernel(mesh->pmlNelements,
                            mesh->o_pmlElementIds,
                            mesh->o_pmlIds,
                            currentTime,
                            intfx, intfy, intfz,
                            mesh->o_sgeo,
                            mesh->o_LIFTT,
                            mesh->o_vmapM,
                            mesh->o_vmapP,
                            mesh->o_EToB,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            bns->o_rkq,
                            bns->o_rhsq,
                            bns->o_pmlrhsqx,
                            bns->o_pmlrhsqy,
                            bns->o_pmlrhsqz);
      occaTimerToc(mesh->device,"PmlSurfaceKernel");
    }

    if(mesh->nonPmlNelements){
      occaTimerTic(mesh->device,"NonPmlSurfaceKernel");
      bns->surfaceKernel(mesh->nonPmlNelements,
                         mesh->o_nonPmlElementIds,
                         currentTime,
                         intfx, intfy, intfz,
                         mesh->o_sgeo,
                         mesh->o_LIFTT,
                         mesh->o_vmapM,
                         mesh->o_vmapP,
                         mesh->o_EToB,
                         mesh->o_x,
                         mesh->o_y,
                         mesh->o_z,
                         bns->o_rkq,
                         bns->o_rhsq);
      occaTimerToc(mesh->device,"NonPmlSurfaceKernel");
    }
    occaTimerToc(mesh->device,"SurfaceKernel");

    
    //UPDATE
    occaTimerTic(mesh->device,"UpdateKernel");


    //printf("running with %d pml Nelements\n",mesh->pmlNelements);    
    if (mesh->pmlNelements){   
      occaTimerTic(mesh->device,"PmlUpdateKernel");
      bns->pmlUpdateKernel(mesh->pmlNelements,
			   mesh->o_pmlElementIds,
			   mesh->o_pmlIds,
			   offset,
			   pmloffset,
			   rk,
			   bns->dt,
			   bns->o_sarkC,
			   bns->o_rkA, 
			   bns->o_rkE,
			   bns->o_sarkA, 
			   bns->o_sarkE,
			   bns->o_q,
			   bns->o_pmlqx,
			   bns->o_pmlqy,
			   bns->o_pmlqz,
			   bns->o_rhsq,
			   bns->o_pmlrhsqx,
			   bns->o_pmlrhsqy,
			   bns->o_pmlrhsqz,
			   bns->o_rkrhsq,
			   bns->o_rkrhsqx,
			   bns->o_rkrhsqy,
			   bns->o_rkrhsqz,
			   bns->o_rkq,
			   bns->o_rkqx,
			   bns->o_rkqy,
			   bns->o_rkqz,
			   bns->o_rkerr);
      occaTimerToc(mesh->device,"PmlUpdateKernel");

    }

    if(mesh->nonPmlNelements){
      occaTimerTic(mesh->device,"NonPmlUpdateKernel");
      bns->updateKernel(mesh->nonPmlNelements,
			mesh->o_nonPmlElementIds,
			offset,
			rk,
			bns->dt,
			bns->o_sarkC,
			bns->o_rkA, 
			bns->o_rkE,
			bns->o_sarkA, 
			bns->o_sarkE,
			bns->o_q,
			bns->o_rhsq,
			bns->o_rkrhsq,
			bns->o_rkq,
			bns->o_rkerr);
      occaTimerToc(mesh->device,"NonPmlUpdateKernel");
    }

    occaTimerToc(mesh->device,"UpdateKernel");
    
  }


}
