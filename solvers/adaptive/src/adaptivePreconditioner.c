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

#include "adaptive.h"

void adaptivePreconditioner(adaptive_t *adaptive, dfloat lambda,
                            occa::memory &o_r, occa::memory &o_z){

  mesh_t *mesh = adaptive->mesh;
  precon_t *precon = adaptive->precon;
  setupAide options = adaptive->options;

  if(options.compareArgs("PRECONDITIONER", "JACOBI")){
    
    dlong Ntotal = mesh->Np*mesh->Nelements;
    // Jacobi preconditioner
    
    adaptive->dotMultiplyKernel(Ntotal, o_r, precon->o_invDiagA, o_z);
  }
  else if (options.compareArgs("PRECONDITIONER", "MULTIGRID")) {

    parAlmond::Precon(precon->parAlmond, o_z, o_r);

  } else if (options.compareArgs("PRECONDITIONER", "FULLALMOND")) {

    if (options.compareArgs("DISCRETIZATION", "IPDG")) {
      parAlmond::Precon(precon->parAlmond, o_z, o_r);
    } else if (options.compareArgs("DISCRETIZATION", "CONTINUOUS")) {
      ogsGather(precon->o_rhsG, o_r, ogsDfloat, ogsAdd, adaptive->ogs);
      adaptive->dotMultiplyKernel(adaptive->ogs->Ngather,
                      adaptive->ogs->o_gatherInvDegree, precon->o_rhsG, precon->o_rhsG);

      parAlmond::Precon(precon->parAlmond, precon->o_xG, precon->o_rhsG);
      ogsScatter(o_z, precon->o_xG, ogsDfloat, ogsAdd, adaptive->ogs);
    }

  } else if(options.compareArgs("PRECONDITIONER", "MASSMATRIX")){

    dfloat invLambda = 1./lambda;

    if (options.compareArgs("DISCRETIZATION", "IPDG")) {
      precon->blockJacobiKernel(mesh->Nelements, invLambda, mesh->o_vgeo, precon->o_invMM, o_r, o_z);
    } else if (options.compareArgs("DISCRETIZATION", "CONTINUOUS")) {
      ogs_t *ogs = adaptive->ogs;

      adaptive->dotMultiplyKernel(mesh->Nelements*mesh->Np, ogs->o_invDegree, o_r, adaptive->o_rtmp);

      if(mesh->NglobalGatherElements)
        precon->partialblockJacobiKernel(mesh->NglobalGatherElements,
                                mesh->o_globalGatherElementList,
                                invLambda, mesh->o_vgeo, precon->o_invMM, adaptive->o_rtmp, o_z);

      ogsGatherScatterStart(o_z, ogsDfloat, ogsAdd, ogs);

      if(mesh->NlocalGatherElements)
        precon->partialblockJacobiKernel(mesh->NlocalGatherElements,
                                mesh->o_localGatherElementList,
                                invLambda, mesh->o_vgeo, precon->o_invMM, adaptive->o_rtmp, o_z);

      ogsGatherScatterFinish(o_z, ogsDfloat, ogsAdd, ogs);

      adaptive->dotMultiplyKernel(mesh->Nelements*mesh->Np, ogs->o_invDegree, o_z, o_z);

      //post-mask
      if (adaptive->Nmasked) mesh->maskKernel(adaptive->Nmasked, adaptive->o_maskIds, o_z);
    }

  } else if (options.compareArgs("PRECONDITIONER", "SEMFEM")) {

    if (adaptive->elementType==TRIANGLES||adaptive->elementType==TETRAHEDRA) {
      o_z.copyFrom(o_r);
      adaptive->dotMultiplyKernel(mesh->Nelements*mesh->Np, adaptive->o_invDegree, o_z, o_z);
      precon->SEMFEMInterpKernel(mesh->Nelements,mesh->o_SEMFEMAnterp,o_z,precon->o_rFEM);
      ogsGather(precon->o_GrFEM, precon->o_rFEM, ogsDfloat, ogsAdd, precon->FEMogs);

      parAlmond::Precon(precon->parAlmond, precon->o_GzFEM, precon->o_GrFEM);

      ogsScatter(precon->o_zFEM, precon->o_GzFEM, ogsDfloat, ogsAdd, precon->FEMogs);
      precon->SEMFEMAnterpKernel(mesh->Nelements,mesh->o_SEMFEMAnterp,precon->o_zFEM,o_z);
      adaptive->dotMultiplyKernel(mesh->Nelements*mesh->Np, adaptive->o_invDegree, o_z, o_z);

      ogsGatherScatter(o_z, ogsDfloat, ogsAdd, adaptive->ogs);
      if (adaptive->Nmasked) mesh->maskKernel(adaptive->Nmasked, adaptive->o_maskIds, o_z);
    } else {
      ogsGather(precon->o_rhsG, o_r, ogsDfloat, ogsAdd, precon->FEMogs);
      adaptive->dotMultiplyKernel(precon->FEMogs->Ngather,
                      precon->FEMogs->o_gatherInvDegree, precon->o_rhsG, precon->o_rhsG);

      parAlmond::Precon(precon->parAlmond, precon->o_xG, precon->o_rhsG);

      ogsScatter(o_z, precon->o_xG, ogsDfloat, ogsAdd, precon->FEMogs);
    }
  }
  else{ // turn off preconditioner
    o_z.copyFrom(o_r);
  }

#if USE_NULL_PROJECTION==1
  if(adaptive->allNeumann) // zero mean of RHS
    adaptiveZeroMean(adaptive, o_z);
#endif
}


