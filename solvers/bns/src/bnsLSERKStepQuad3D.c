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

    }

    // COMPUTE RAMP FUNCTION 
    dfloat fx=0, fy=0, fz=0, intfx=0, intfy=0, intfz=0;
    //    bnsBodyForce(t, &fx, &fy, &fz, &intfx, &intfy, &intfz);

    // compute volume contribution to DG boltzmann RHS added d/dt (ramp(qbar)) to RHS
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

#if 1
    // compute relaxation terms using cubature
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
#endif
    
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);
      // copy halo data to DEVICE
      size_t offset = mesh->Np*bns->Nfields*mesh->Nelements*sizeof(dfloat); // offset for halo data
      bns->o_q.copyFrom(recvBuffer, haloBytes, offset);
    }



    // SURFACE KERNELS
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
    
    //UPDATE
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
}
