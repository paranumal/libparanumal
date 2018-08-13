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

// compute LU = L(U)
void insDiffusion(ins_t *ins, dfloat time, occa::memory o_U, occa::memory o_LU){

  mesh_t *mesh = ins->mesh;
  elliptic_t *uSolver = ins->uSolver; //borrow the uSolver for the gather lists
  setupAide options = ins->vOptions;

  if(options.compareArgs("DISCRETIZATION", "CONTINUOUS")){
    ogs_t *ogs = mesh->ogs;

    if(mesh->NglobalGatherElements)
      ins->diffusionKernel(mesh->NglobalGatherElements, 
                           mesh->o_globalGatherElementList,
                           mesh->o_ggeo, 
                           mesh->o_vgeo, 
                           mesh->o_sgeo, 
                           mesh->o_Dmatrices, 
                           mesh->o_Smatrices, 
                           mesh->o_vmapM,
                           mesh->o_sMT,
                           ins->nu,
                           time,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_z,
                           ins->o_VmapB,
                           ins->fieldOffset,
                           o_U, 
                           o_LU);

    ogsGatherScatterManyStart(o_LU, ins->dim, ins->fieldOffset,
                              ogsDfloat, ogsAdd, ogs);

    if(mesh->NlocalGatherElements)
        ins->diffusionKernel(mesh->NlocalGatherElements, 
                             mesh->o_localGatherElementList,
                             mesh->o_ggeo, 
                             mesh->o_vgeo, 
                             mesh->o_sgeo, 
                             mesh->o_Dmatrices, 
                             mesh->o_Smatrices, 
                             mesh->o_vmapM,
                             mesh->o_sMT,
                             ins->nu,
                             time,
                             mesh->o_x,
                             mesh->o_y,
                             mesh->o_z,
                             ins->o_VmapB,
                             ins->fieldOffset,
                             o_U, 
                             o_LU);

    ogsGatherScatterManyFinish(o_LU, ins->dim, ins->fieldOffset,
                              ogsDfloat, ogsAdd, ogs);    

  } else if(options.compareArgs("DISCRETIZATION", "IPDG")) {
    dlong offset = 0;

    if(mesh->totalHaloPairs>0){
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
    
    ins->velocityGradientKernel(mesh->Nelements,
                                offset,
                                mesh->o_vgeo,
                                mesh->o_Dmatrices,
                                ins->fieldOffset,
                                o_U,
                                ins->o_GU);

    if(mesh->NinternalElements) {
      ins->diffusionIpdgKernel(mesh->NinternalElements,
                               mesh->o_internalElementIds,
                               mesh->o_vmapM,
                               mesh->o_vmapP,
                               ins->nu,
                               uSolver->tau,
                               mesh->o_vgeo,
                               mesh->o_sgeo,
                               mesh->o_EToB,
                               time,
                               mesh->o_x,
                               mesh->o_y,
                               mesh->o_z,
                               ins->fieldOffset,
                               mesh->o_Dmatrices,
                               mesh->o_LIFTT,
                               ins->o_GU,
                               o_LU);
    }

    // COMPLETE HALO EXCHANGE
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);

      ins->o_vHaloBuffer.copyFrom(ins->vRecvBuffer); 

      ins->velocityHaloScatterKernel(mesh->Nelements,
                                    mesh->totalHaloPairs,
                                    ins->fieldOffset,
                                    o_U,
                                    ins->o_vHaloBuffer);
    }

    if(mesh->totalHaloPairs){
      offset = mesh->Nelements;  
      ins->velocityGradientKernel(mesh->totalHaloPairs,
                                  offset,
                                  mesh->o_vgeo,
                                  mesh->o_Dmatrices,
                                  ins->fieldOffset,
                                  o_U,
                                  ins->o_GU);
    }

    if(mesh->NnotInternalElements) {
      ins->diffusionIpdgKernel(mesh->NnotInternalElements,
                              mesh->o_notInternalElementIds,
                              mesh->o_vmapM,
                              mesh->o_vmapP,
                              ins->nu,
                              uSolver->tau,
                              mesh->o_vgeo,
                              mesh->o_sgeo,
                              mesh->o_EToB,
                              time,
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              ins->fieldOffset,
                              mesh->o_Dmatrices,
                              mesh->o_LIFTT,
                              ins->o_GU,
                              o_LU);
    }
  } 
}