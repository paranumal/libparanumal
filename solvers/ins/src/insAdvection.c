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

#include "ins.h"

// compute NU = N(U)
void insAdvection(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_NU){

  mesh_t *mesh = ins->mesh;
  
  //Exctract Halo On Device, all fields
  if(mesh->totalHaloPairs>0 && !ins->solveHeat){
    ins->velocityHaloExtractKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 mesh->o_haloElementList,
                                 ins->fieldOffset,
                                 o_U,
                                 ins->o_vHaloBuffer);

    // copy extracted halo to HOST 
    ins->o_vHaloBuffer.copyTo(ins->vSendBuffer);           
  
    // start halo exchange
    meshHaloExchangeStart(mesh,
                         mesh->Np*(ins->NVfields)*sizeof(dfloat),
                         ins->vSendBuffer,
                         ins->vRecvBuffer);
  }

  // Compute Volume Contribution
  occaTimerTic(mesh->device,"AdvectionVolume");
  if(ins->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
    ins->advectionCubatureVolumeKernel(mesh->Nelements,
                                       mesh->o_vgeo,
                                       mesh->o_cubvgeo,
                                       mesh->o_cubDWmatrices,
                                       mesh->o_cubInterpT,
                                       mesh->o_cubProjectT,
                                       ins->fieldOffset,
                                       o_U,
                                       ins->o_cU,
                                       o_NU);
  } else {
    ins->advectionVolumeKernel(mesh->Nelements,
                               mesh->o_vgeo,
                               mesh->o_Dmatrices,
                               ins->fieldOffset,
                               o_U,
                               o_NU);
  }
  occaTimerToc(mesh->device,"AdvectionVolume");

  // COMPLETE HALO EXCHANGE
  if(mesh->totalHaloPairs>0 && !ins->solveHeat){
    meshHaloExchangeFinish(mesh);

    ins->o_vHaloBuffer.copyFrom(ins->vRecvBuffer); 

    ins->velocityHaloScatterKernel(mesh->Nelements,
                                  mesh->totalHaloPairs,
                                  ins->fieldOffset,
                                  o_U,
                                  ins->o_vHaloBuffer);
  }

  occaTimerTic(mesh->device,"AdvectionSurface");
  if(ins->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
    ins->advectionCubatureSurfaceKernel(mesh->Nelements,
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
                                        ins->fieldOffset,
                                        o_U,
                                        o_NU);
  } else {
    ins->advectionSurfaceKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_LIFTT,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                time,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                ins->fieldOffset,
                                o_U,
                                o_NU);
  }
  occaTimerToc(mesh->device,"AdvectionSurface");
}
