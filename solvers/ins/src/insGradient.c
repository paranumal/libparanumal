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

// complete a time step using LSERK4
void insGradient(ins_t *ins, dfloat time, occa::memory o_P, occa::memory o_GP){

  mesh_t *mesh = ins->mesh;
  
  if (ins->pOptions.compareArgs("DISCRETIZATION","IPDG")) {
    if(mesh->totalHaloPairs>0){
      ins->pressureHaloExtractKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 mesh->o_haloElementList,
                                 o_P,
                                 ins->o_pHaloBuffer);

      // copy extracted halo to HOST
      ins->o_pHaloBuffer.copyTo(ins->pSendBuffer);

      // start halo exchange
      meshHaloExchangeStart(mesh,
                           mesh->Np*sizeof(dfloat),
                           ins->pSendBuffer,
                           ins->pRecvBuffer);
    }
  }

  occaTimerTic(mesh->device,"GradientVolume");
  // Compute Volume Contribution
  ins->gradientVolumeKernel(mesh->Nelements,
                            mesh->o_vgeo,
                            mesh->o_Dmatrices,
                            ins->fieldOffset,
                            o_P,
                            o_GP);
  occaTimerToc(mesh->device,"GradientVolume");

  // COMPLETE HALO EXCHANGE
  if (ins->pOptions.compareArgs("DISCRETIZATION","IPDG")) {
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);

      ins->o_pHaloBuffer.copyFrom(ins->pRecvBuffer);

      ins->pressureHaloScatterKernel(mesh->Nelements,
                                      mesh->totalHaloPairs,
                                      o_P,
                                      ins->o_pHaloBuffer);
    }

    occaTimerTic(mesh->device,"GradientSurface");
    // Compute Surface Conribution
    ins->gradientSurfaceKernel(mesh->Nelements,
                               mesh->o_sgeo,
                               mesh->o_LIFTT,
                               mesh->o_vmapM,
                               mesh->o_vmapP,
                               mesh->o_EToB,
                               mesh->o_x,
                               mesh->o_y,
                               mesh->o_z,
                               time,
                               ins->fieldOffset,
                               o_P,
                               o_GP);
    occaTimerToc(mesh->device,"GradientSurface");
  }
}
