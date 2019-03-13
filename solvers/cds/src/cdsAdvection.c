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
#define ADV_DEBUG 0

#include "cds.h"

// compute NS = N(UxS)
void cdsAdvection(cds_t *cds, dfloat time, occa::memory o_U, occa::memory o_S, occa::memory o_NS){

#if ADV_DEBUG
  printf("Starting Advection Step \n "); 
#endif
  mesh_t *mesh = cds->mesh;

   //Exctract Halo On Device, all fields
  if(mesh->totalHaloPairs>0){
#if 0
    cds->haloExtractKernel(mesh->Nelements,
                           mesh->totalHaloPairs,
                           mesh->o_haloElementList,
                           cds->vOffset,
                           cds->sOffset,
                           o_U,
                           o_S,
                           cds->o_haloBuffer);

    // copy extracted halo to HOST 
    cds->o_haloBuffer.copyTo(cds->sendBuffer);           
  
    // start halo exchange
    meshHaloExchangeStart(mesh,
        mesh->Np*(cds->NSfields+cds->NVfields)*sizeof(dfloat),
        cds->sendBuffer,
        cds->recvBuffer);

#else
      mesh->device.finish();
      // switch to data stream
      mesh->device.setStream(mesh->dataStream);

      cds->haloGetKernel(mesh->totalHaloPairs,
                         cds->vOffset,
                         cds->sOffset,
                         mesh->o_haloElementList,
                         mesh->o_haloGetNodeIds,
                         o_U,
                         o_S,
                         cds->o_haloBuffer);
          
      dlong Ndata = (cds->NVfields+cds->NSfields)*mesh->Nfp*mesh->totalHaloPairs;
          // copy extracted halo to HOST 
      cds->o_haloBuffer.copyTo(cds->sendBuffer, Ndata*sizeof(dfloat), 0, "async: true");// zero offset

      mesh->device.setStream(mesh->defaultStream);

#endif
  }

#if ADV_DEBUG
  printf("\tSolving  Advection Volume ..."); 
#endif

  
  // Compute Volume Contribution
  occaTimerTic(mesh->device,"AdvectionVolume");
  if(cds->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
    cds->advectionCubatureVolumeKernel(mesh->Nelements,
                                       mesh->o_vgeo,
                                       mesh->o_cubvgeo,
                                       mesh->o_cubDWmatrices,
                                       mesh->o_cubInterpT,
                                       mesh->o_cubProjectT,
                                       cds->vOffset,
                                       cds->sOffset,
				                               o_U,
                                       o_S,
                                       o_NS);
  } else {
    cds->advectionVolumeKernel(mesh->Nelements,
                               mesh->o_vgeo,
                               mesh->o_Dmatrices,
                               cds->vOffset,
                               cds->sOffset,
                               o_U,
                               o_S,
                               o_NS);
  }
  occaTimerToc(mesh->device,"AdvectionVolume");

#if ADV_DEBUG
  printf("\tdone \n "); 
#endif

  // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){

#if 0
    meshHaloExchangeFinish(mesh);

    cds->o_haloBuffer.copyFrom(cds->recvBuffer); 

    cds->haloScatterKernel(mesh->Nelements,
			   mesh->totalHaloPairs,
         cds->vOffset,
         cds->sOffset,
         o_U,
         o_S,
			   cds->o_haloBuffer);

     // cds->haloScatterKernel(mesh->Nelements,
     //     mesh->totalHaloPairs,
     //     cds->sOffset,
     //     o_S,
     //     cds->o_haloBuffer);

#else
      // make sure compute device is ready to perform halo extract
      mesh->device.setStream(mesh->dataStream);
      mesh->device.finish();

      // start halo exchange
      meshHaloExchangeStart(mesh,
                            mesh->Nfp*(cds->NVfields+ cds->NSfields)*sizeof(dfloat),
                            cds->sendBuffer,
                            cds->recvBuffer);

      meshHaloExchangeFinish(mesh);

      dlong Ndata = (cds->NVfields + cds->NSfields)*mesh->Nfp*mesh->totalHaloPairs;
      cds->o_haloBuffer.copyFrom(cds->recvBuffer, Ndata*sizeof(dfloat), 0, "async: true");  

      cds->haloPutKernel(mesh->totalHaloPairs,
                          cds->vOffset,
                          cds->sOffset,
                          mesh->o_haloElementList,
                          mesh->o_haloPutNodeIds,
                          o_U,
                          o_S,
                          cds->o_haloBuffer);

      mesh->device.finish();

      mesh->device.setStream(mesh->defaultStream);
      mesh->device.finish();
#endif
  }


#if ADV_DEBUG
  printf("\tSolving  Advection Surface ..."); 
#endif


  
  occaTimerTic(mesh->device,"AdvectionSurface");
  if(cds->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
    cds->advectionCubatureSurfaceKernel(mesh->Nelements,
                                        mesh->o_vgeo,
                                        mesh->o_sgeo,
                                        mesh->o_cubsgeo,
                                        mesh->o_intInterpT,
                                        mesh->o_intLIFTT,
                                        mesh->o_cubInterpT,
                                        mesh->o_cubProjectT,
                                        mesh->o_vmapM,
                                        mesh->o_vmapP,
                                        mesh->o_EToB,
                                        time,
                                        mesh->o_intx,
                                        mesh->o_inty,
                                        mesh->o_intz,
                                        cds->vOffset,
                                        cds->sOffset,
                                        o_U,
                                        o_S,
                                        o_NS);
  } else {
    cds->advectionSurfaceKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_LIFTT,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                time,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                cds->vOffset,
				                        cds->sOffset,
                                o_U,
				                        o_S,
                                o_NS);
  }
  occaTimerToc(mesh->device,"AdvectionSurface");

#if ADV_DEBUG
  printf("\tdone\n"); 
#endif
  
#if ADV_DEBUG
  printf("Advection Step Completed\n"); 
#endif

}
