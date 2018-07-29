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

    if(ogs->NhaloGather) {
      ins->diffusionKernel(uSolver->NglobalGatherElements, 
                           uSolver->o_globalGatherElementList,
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

      mesh->device.finish();
      mesh->device.setStream(uSolver->dataStream);
      mesh->gatherKernel(ogs->NhaloGather, 
                         ogs->o_haloGatherOffsets, 
                         ogs->o_haloGatherLocalIds, 
                         o_LU, 
                         ins->NVfields,
                         ins->fieldOffset,
                         ins->o_velocityHaloGatherTmp);
      ins->o_velocityHaloGatherTmp.copyTo(ins->velocityHaloGatherTmp,"async: true");
      mesh->device.setStream(uSolver->defaultStream);
    }
    if(uSolver->NlocalGatherElements){
        ins->diffusionKernel(uSolver->NlocalGatherElements, 
                             uSolver->o_localGatherElementList,
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
    }

    // finalize gather using local and global contributions
    if(ogs->NnonHaloGather) 
      mesh->gatherScatterKernel(ogs->NnonHaloGather, 
                                ogs->o_nonHaloGatherOffsets, 
                                ogs->o_nonHaloGatherLocalIds,
                                ins->NVfields,
                                ins->fieldOffset, 
                                o_LU);

    // C0 halo gather-scatter (on data stream)
    if(ogs->NhaloGather) {
      mesh->device.setStream(uSolver->dataStream);
      mesh->device.finish();

      // MPI based gather scatter using libgs
      gsVecParallelGatherScatter(ogs->haloGsh, ins->velocityHaloGatherTmp, ins->NVfields, dfloatString, "add");

      // copy totally gather halo data back from HOST to DEVICE
      ins->o_velocityHaloGatherTmp.copyFrom(ins->velocityHaloGatherTmp,"async: true");
    
      // do scatter back to local nodes
      mesh->scatterKernel(ogs->NhaloGather, 
                          ogs->o_haloGatherOffsets, 
                          ogs->o_haloGatherLocalIds, 
                          ins->NVfields,
                          ins->fieldOffset,
                          ins->o_velocityHaloGatherTmp, 
                          o_LU);
      mesh->device.setStream(uSolver->defaultStream);
    }

    mesh->device.finish();    
    mesh->device.setStream(uSolver->dataStream);
    mesh->device.finish();    
    mesh->device.setStream(uSolver->defaultStream);

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