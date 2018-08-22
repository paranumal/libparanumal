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

#include "mppf.h"

// solve lambda*U + A*U = rhsU
void mppfVelocitySolve(mppf_t *mppf, dfloat time, occa::memory o_rkU){
  
  mesh_t *mesh = mppf->mesh; 
  elliptic_t *usolver = mppf->uSolver; 
  elliptic_t *vsolver = mppf->vSolver; 
  elliptic_t *wsolver = mppf->wSolver; 
  
  if (mppf->vOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    mppf->velocityRhsBCKernel(mesh->Nelements,
                              mesh->o_ggeo,
                              mesh->o_sgeo,
                              mesh->o_Dmatrices,
                              mesh->o_Smatrices,
                              mesh->o_MM,
                              mesh->o_vmapM,
			                        mesh->o_EToB,
                              mesh->o_sMT,
                              mppf->lambdaVel,
                              time,
                              mesh->o_x,
                              mesh->o_y,
                              mesh->o_z,
                              mppf->o_VmapB,
                              mppf->o_rhsU,
                              mppf->o_rhsV,
                              mppf->o_rhsW);
    
    // gather-scatter
    ellipticParallelGatherScatter(mesh, mesh->ogs, mppf->o_rhsU, dfloatString, "add");  
    ellipticParallelGatherScatter(mesh, mesh->ogs, mppf->o_rhsV, dfloatString, "add");  
    if (mppf->dim==3)
      ellipticParallelGatherScatter(mesh, mesh->ogs, mppf->o_rhsW, dfloatString, "add");  
    
    if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, mppf->o_rhsU);
    if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, mppf->o_rhsV);
    if (mppf->dim==3)
      if (wsolver->Nmasked) mesh->maskKernel(wsolver->Nmasked, wsolver->o_maskIds, mppf->o_rhsW);

  } else if (mppf->vOptions.compareArgs("DISCRETIZATION","IPDG")) {

    occaTimerTic(mesh->device,"velocityRhsIpdg");    
    mppf->velocityRhsIpdgBCKernel(mesh->Nelements,
                                  mesh->o_vmapM,
                                  usolver->tau,
                                  time,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_Dmatrices,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  mppf->o_rhsU,
                                  mppf->o_rhsV,
                                  mppf->o_rhsW);
    occaTimerToc(mesh->device,"velocityRhsIpdg");   
  }

  //copy current velocity fields as initial guess? (could use Uhat or beter guess)
  dlong Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  mppf->o_UH.copyFrom(mppf->o_U,Ntotal*sizeof(dfloat),0,0*mppf->fieldOffset*sizeof(dfloat));
  mppf->o_VH.copyFrom(mppf->o_U,Ntotal*sizeof(dfloat),0,1*mppf->fieldOffset*sizeof(dfloat));
  if (mppf->dim==3)
    mppf->o_WH.copyFrom(mppf->o_U,Ntotal*sizeof(dfloat),0,2*mppf->fieldOffset*sizeof(dfloat));

  if (mppf->vOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, mppf->o_UH);
    if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, mppf->o_VH);
    if (mppf->dim==3)
      if (wsolver->Nmasked) mesh->maskKernel(wsolver->Nmasked, wsolver->o_maskIds, mppf->o_WH);

  }
  
  occaTimerTic(mesh->device,"Ux-Solve");
  mppf->NiterU = ellipticSolve(usolver, mppf->lambdaVel, mppf->velTOL, mppf->o_rhsU, mppf->o_UH);
  occaTimerToc(mesh->device,"Ux-Solve"); 

  occaTimerTic(mesh->device,"Uy-Solve");
  mppf->NiterV = ellipticSolve(vsolver, mppf->lambdaVel, mppf->velTOL, mppf->o_rhsV, mppf->o_VH);
  occaTimerToc(mesh->device,"Uy-Solve");

  if (mppf->dim==3) {
    occaTimerTic(mesh->device,"Uz-Solve");
    mppf->NiterW = ellipticSolve(wsolver, mppf->lambdaVel, mppf->velTOL, mppf->o_rhsW, mppf->o_WH);
    occaTimerToc(mesh->device,"Uz-Solve");
  }

  if (mppf->vOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    mppf->velocityAddBCKernel(mesh->Nelements,
                            time,
                            mesh->o_sgeo,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            mppf->o_VmapB,
                            mppf->o_UH,
                            mppf->o_VH,
                            mppf->o_WH);
  }

  //copy into intermediate stage storage
  mppf->o_UH.copyTo(o_rkU,Ntotal*sizeof(dfloat),0*mppf->fieldOffset*sizeof(dfloat),0);
  mppf->o_VH.copyTo(o_rkU,Ntotal*sizeof(dfloat),1*mppf->fieldOffset*sizeof(dfloat),0);    
  if (mppf->dim==3)
    mppf->o_WH.copyTo(o_rkU,Ntotal*sizeof(dfloat),2*mppf->fieldOffset*sizeof(dfloat),0);    
}
