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

#include "cds.h"

// compute NS = N(UxS)
void cdsAdvection(cds_t *cds, dfloat time, occa::memory o_U, occa::memory o_S, occa::memory o_NS){

  mesh_t *mesh = cds->mesh;
  
  //Exctract Halo On Device, all fields
  if(mesh->totalHaloPairs>0){
    cds->haloExtractKernel(mesh->Nelements,
                           mesh->totalHaloPairs,
                           mesh->o_haloElementList,
                           cds->fieldOffset,
                           cds->o_S,
                           cds->o_haloBuffer);

    // copy extracted halo to HOST 
    cds->o_haloBuffer.copyTo(cds->sendBuffer);           
  
    // start halo exchange
    meshHaloExchangeStart(mesh,
                         mesh->Np*(cds->NSfields)*sizeof(dfloat),
                         cds->sendBuffer,
                         cds->recvBuffer);
  }


  /*
  // Compute Volume Contribution
  occaTimerTic(mesh->device,"AdvectionVolume");
  if(cds->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
    cds->advectionCubatureVolumeKernel(mesh->Nelements,
                                       mesh->o_vgeo,
                                       mesh->o_cubvgeo,
                                       mesh->o_cubDWmatrices,
                                       mesh->o_cubInterpT,
                                       mesh->o_cubProjectT,
                                       cds->fieldOffset,
                                       o_U,
                                       cds->o_cU,
                                       o_NU);
  } else {
    cds->advectionVolumeKernel(mesh->Nelements,
                               mesh->o_vgeo,
                               mesh->o_Dmatrices,
                               cds->fieldOffset,
                               o_U,
                               o_NU);
  }
  occaTimerToc(mesh->device,"AdvectionVolume");

  // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0){
    meshHaloExchangeFinish(mesh);

    cds->o_vHaloBuffer.copyFrom(cds->vRecvBuffer); 

    cds->velocityHaloScatterKernel(mesh->Nelements,
                                  mesh->totalHaloPairs,
                                  cds->fieldOffset,
                                  o_U,
                                  cds->o_vHaloBuffer);
  }

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
                                        cds->fieldOffset,
                                        o_U,
                                        o_NU);
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
                                cds->fieldOffset,
                                o_U,
                                o_NU);
  }
  occaTimerToc(mesh->device,"AdvectionSurface");

  */
}
