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
void bnsLSERKStep(bns_t *bns, int tstep, int haloBytes,
		  dfloat * sendBuffer, dfloat *recvBuffer, setupAide &options){


  const dlong offset    = 0.0;
  const dlong pmloffset = 0.0;
  const int shift       = 0; 
  // bns->shiftIndex = 0; 
  mesh_t *mesh = bns->mesh;   

  // LSERK4 stages
  for(int rk=0;rk<mesh->Nrk;++rk){

    // intermediate stage time
    dfloat t = bns->startTime + tstep*bns->dt + bns->dt*mesh->rkc[rk];

    if(mesh->totalHaloPairs>0){
#if BNS_ASYNC 
      mesh->device.setStream(dataStream);

      int Nentries = mesh->Np*bns->Nfields;
      mesh->haloExtractKernel(mesh->totalHaloPairs,
			      Nentries,
			      mesh->o_haloElementList,
			      bns->o_q,
			      mesh->o_haloBuffer);

      // copy extracted halo to HOST
      mesh->o_haloBuffer.copyTo(sendBuffer,"async: true");
      mesh->device.setStream(defaultStream);

#else
      int Nentries = mesh->Np*bns->Nfields;
      mesh->haloExtractKernel(mesh->totalHaloPairs,
			      Nentries,
			      mesh->o_haloElementList,
			      bns->o_q,
			      mesh->o_haloBuffer);
      // copy extracted halo to HOST
      mesh->o_haloBuffer.copyTo(sendBuffer);
      // start halo exchange
      meshHaloExchangeStart(mesh, bns->Nfields*mesh->Np*sizeof(dfloat), sendBuffer, recvBuffer);
#endif
    }

    // COMPUTE RAMP FUNCTION 
    dfloat fx, fy, fz, intfx, intfy, intfz;
    bnsBodyForce(t, &fx, &fy, &fz, &intfx, &intfy, &intfz);

    occaTimerTic(mesh->device, "VolumeKernel");    
    // compute volume contribution to DG boltzmann RHS
    if(mesh->pmlNelements){ 
      occaTimerTic(mesh->device,"PmlVolumeKernel");

      if(bns->pmlcubature){
        bns->pmlVolumeKernel(mesh->pmlNelements,
			     mesh->o_pmlElementIds,
			     mesh->o_pmlIds,
			     offset,    // 0
			     pmloffset, // 0
			     shift, // 0
			     fx, fy, fz, 
			     mesh->o_vgeo,
			     mesh->o_Dmatrices,
			     bns->o_q,
			     bns->o_pmlqx,
			     bns->o_pmlqy,
			     bns->o_pmlqz,
			     bns->o_rhsq,
			     bns->o_pmlrhsqx,
			     bns->o_pmlrhsqy,
			     bns->o_pmlrhsqz);
      }else{
	bns->pmlVolumeKernel(mesh->pmlNelements,
			     mesh->o_pmlElementIds,
			     mesh->o_pmlIds,
			     offset,    // 0
			     pmloffset, // 0
			     shift, 
			     fx,fy,fz,
			     mesh->o_vgeo,
			     bns->o_pmlSigmaX,
			     bns->o_pmlSigmaY,
			     bns->o_pmlSigmaZ,
			     mesh->o_Dmatrices,
			     bns->o_q,
			     bns->o_pmlqx,
			     bns->o_pmlqy,
			     bns->o_pmlqz,
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
			offset, 
			shift,
			fx, fy, fz,
			mesh->o_vgeo,
			mesh->o_x,
			mesh->o_y,
			mesh->o_z,
			mesh->o_Dmatrices,
			bns->o_q,
			bns->o_rhsq);
      occaTimerToc(mesh->device,"NonPmlVolumeKernel");
    }
    occaTimerToc(mesh->device, "VolumeKernel");    
    


#if 1
    // TW: TURN OFF RELAXATION FOR TESTING ONLY
    occaTimerTic(mesh->device, "RelaxationKernel");
    if(mesh->pmlNelements){
      occaTimerTic(mesh->device, "PmlRelaxationKernel");

      if(bns->pmlcubature){
	bns->pmlRelaxationKernel(mesh->pmlNelements,
				 mesh->o_pmlElementIds,
				 mesh->o_pmlIds,
				 mesh->o_vgeo,
				 mesh->o_cubvgeo,
				 offset,    // 0
				 pmloffset, // 0
				 shift, // 0
				 mesh->o_cubInterpT,
				 mesh->o_cubProjectT,
				 bns->o_pmlSigmaX,
				 bns->o_pmlSigmaY,
				 bns->o_pmlSigmaZ,
				 bns->o_q,
				 bns->o_pmlqx,
				 bns->o_pmlqy,
				 bns->o_pmlqz,
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
				 offset,    // 0
				 pmloffset, // 0
				 shift, // 0
				 mesh->o_cubInterpT,
				 mesh->o_cubProjectT,
				 bns->o_q,
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
			    offset, // 0
			    shift, //0
			    mesh->o_cubInterpT,
			    mesh->o_cubProjectT,
			    bns->o_q,
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
      bns->o_q.copyFrom(recvBuffer, haloBytes, offset,"async: true");
      mesh->device.finish();

      mesh->device.setStream(defaultStream);
#else
      meshHaloExchangeFinish(mesh);
      // copy halo data to DEVICE
      size_t offset = mesh->Np*bns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      bns->o_q.copyFrom(recvBuffer, haloBytes, offset);
#endif
    }



    // SURFACE KERNELS
    occaTimerTic(mesh->device,"SurfaceKernel");

    if(mesh->pmlNelements){
      occaTimerTic(mesh->device,"PmlSurfaceKernel");
      bns->pmlSurfaceKernel(mesh->pmlNelements,
                            mesh->o_pmlElementIds,
                            mesh->o_pmlIds,
                            t,
                            intfx, intfy, intfz,
                            mesh->o_sgeo,
                            mesh->o_LIFTT,
                            mesh->o_vmapM,
                            mesh->o_vmapP,
                            mesh->o_EToB,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            bns->o_q,
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
			 t,
			 intfx, intfy, intfz,
			 mesh->o_sgeo,
			 mesh->o_LIFTT,
			 mesh->o_vmapM,
			 mesh->o_vmapP,
			 mesh->o_EToB,
			 mesh->o_x,
			 mesh->o_y,
			 mesh->o_z,
			 bns->o_q,
			 bns->o_rhsq);
      occaTimerToc(mesh->device,"NonPmlSurfaceKernel");
    }
    occaTimerToc(mesh->device,"SurfaceKernel");



    if(bns->elementType==QUADRILATERALS && mesh->dim==3){
      bns->constrainKernel(mesh->Nelements, mesh->o_x, mesh->o_y, mesh->o_z, bns->o_rhsq);
    }
    
    
    // ramp function for flow at next RK stage
    dfloat tupdate = tstep*bns->dt + bns->dt*mesh->rkc[rk+1];

    dfloat fxUpdate, fyUpdate, fzUpdate, intfxUpdate, intfyUpdate, intfzUpdate;
    bnsBodyForce(tupdate, &fxUpdate, &fyUpdate, &fzUpdate, &intfxUpdate, &intfyUpdate, &intfzUpdate);
    
    //UPDATE
    occaTimerTic(mesh->device,"UpdateKernel");

    if (mesh->pmlNelements){   
      occaTimerTic(mesh->device,"PmlUpdateKernel");
      bns->pmlUpdateKernel(mesh->pmlNelements,
			   mesh->o_pmlElementIds,
			   mesh->o_pmlIds,
			   bns->dt,
			   mesh->rka[rk],
			   mesh->rkb[rk],
			   bns->o_rhsq,
			   bns->o_pmlrhsqx,
			   bns->o_pmlrhsqy,
			   bns->o_pmlrhsqz,
			   bns->o_resq,
			   bns->o_pmlresqx,
			   bns->o_pmlresqy,
			   bns->o_pmlresqz,
			   bns->o_pmlqx,
			   bns->o_pmlqy,
			   bns->o_pmlqz,
			   bns->o_q);
      occaTimerToc(mesh->device,"PmlUpdateKernel");

    }

    if(mesh->nonPmlNelements){
      occaTimerTic(mesh->device,"NonPmlUpdateKernel");
      bns->updateKernel(mesh->nonPmlNelements,
                        mesh->o_nonPmlElementIds,
                        bns->dt,
                        mesh->rka[rk],
                        mesh->rkb[rk],
                        bns->o_rhsq,
                        bns->o_resq,
                        bns->o_q);
      occaTimerToc(mesh->device,"NonPmlUpdateKernel");
    }

    occaTimerToc(mesh->device,"UpdateKernel");
    
  }
}
