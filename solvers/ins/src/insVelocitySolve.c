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

// solve lambda*U + A*U = rhsU
void insVelocitySolve(ins_t *ins, dfloat time, int stage,  occa::memory o_rhsU, 
                                                           occa::memory o_rhsV, 
                                                           occa::memory o_rhsW, 
                                                           occa::memory o_Uhat){
  
  mesh_t *mesh = ins->mesh; 
  elliptic_t *usolver = ins->uSolver; 
  elliptic_t *vsolver = ins->vSolver; 
  elliptic_t *wsolver = ins->wSolver;

  int quad3D = (ins->dim==3 && ins->elementType==QUADRILATERALS) ? 1 : 0;  
  
  if (ins->vOptions.compareArgs("DISCRETIZATION","CONTINUOUS")){

    if(!quad3D) 
      ins->velocityRhsBCKernel(mesh->Nelements,
                                mesh->o_ggeo,
                                mesh->o_sgeo,
                                mesh->o_Dmatrices,
                                mesh->o_Smatrices,
                                mesh->o_MM,
                                mesh->o_vmapM,
    	                          mesh->o_EToB,
                                mesh->o_sMT,
                                ins->lambda,
                                time,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                ins->o_VmapB,
                                o_rhsU,
                                o_rhsV,
                                o_rhsW);

    // gather-scatter
    ogsGatherScatter(o_rhsU, ogsDfloat, ogsAdd, mesh->ogs);
    ogsGatherScatter(o_rhsV, ogsDfloat, ogsAdd, mesh->ogs);

    if (ins->dim==3)
      ogsGatherScatter(o_rhsW, ogsDfloat, ogsAdd, mesh->ogs);

    if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, o_rhsU);
    if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, o_rhsV);
    if (ins->dim==3)
      if (wsolver->Nmasked) mesh->maskKernel(wsolver->Nmasked, wsolver->o_maskIds, o_rhsW);
    
  } else if (ins->vOptions.compareArgs("DISCRETIZATION","IPDG") && !quad3D) {

    occaTimerTic(mesh->device,"velocityRhsIpdg");    
    ins->velocityRhsIpdgBCKernel(mesh->Nelements,
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
                                  o_rhsU,
                                  o_rhsV,
                                  o_rhsW);
    occaTimerToc(mesh->device,"velocityRhsIpdg");   
  }

  //copy current velocity fields as initial guess? (could use Uhat or beter guess)
  dlong Ntotal = (mesh->Nelements+mesh->totalHaloPairs)*mesh->Np;
  ins->o_UH.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,0*ins->fieldOffset*sizeof(dfloat));
  ins->o_VH.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,1*ins->fieldOffset*sizeof(dfloat));
  if (ins->dim==3)
    ins->o_WH.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,2*ins->fieldOffset*sizeof(dfloat));

  if (ins->vOptions.compareArgs("DISCRETIZATION","CONTINUOUS") && !quad3D) {
    if (usolver->Nmasked) mesh->maskKernel(usolver->Nmasked, usolver->o_maskIds, ins->o_UH);
    if (vsolver->Nmasked) mesh->maskKernel(vsolver->Nmasked, vsolver->o_maskIds, ins->o_VH);
    if (ins->dim==3)
      if (wsolver->Nmasked) mesh->maskKernel(wsolver->Nmasked, wsolver->o_maskIds, ins->o_WH);

  }
  
  occaTimerTic(mesh->device,"Ux-Solve");
  ins->NiterU = ellipticSolve(usolver, ins->lambda, ins->velTOL, o_rhsU, ins->o_UH);
  occaTimerToc(mesh->device,"Ux-Solve"); 

  occaTimerTic(mesh->device,"Uy-Solve");
  ins->NiterV = ellipticSolve(vsolver, ins->lambda, ins->velTOL, o_rhsV, ins->o_VH);
  occaTimerToc(mesh->device,"Uy-Solve");

  if (ins->dim==3) {
    occaTimerTic(mesh->device,"Uz-Solve");
    ins->NiterW = ellipticSolve(wsolver, ins->lambda, ins->velTOL, o_rhsW, ins->o_WH);
    occaTimerToc(mesh->device,"Uz-Solve");
  }

  if (ins->vOptions.compareArgs("DISCRETIZATION","CONTINUOUS") && !quad3D) {
    ins->velocityAddBCKernel(mesh->Nelements,
                            time,
                            mesh->o_sgeo,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            ins->o_VmapB,
                            ins->o_UH,
                            ins->o_VH,
                            ins->o_WH);
  }

  //copy into intermediate stage storage
  ins->o_UH.copyTo(o_Uhat,Ntotal*sizeof(dfloat),0*ins->fieldOffset*sizeof(dfloat),0);
  ins->o_VH.copyTo(o_Uhat,Ntotal*sizeof(dfloat),1*ins->fieldOffset*sizeof(dfloat),0);    
  if (ins->dim==3)
    ins->o_WH.copyTo(o_Uhat,Ntotal*sizeof(dfloat),2*ins->fieldOffset*sizeof(dfloat),0);    
}
