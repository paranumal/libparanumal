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
void insHeatSolve(ins_t *ins, dfloat time, int stage){
  
  mesh_t *mesh = ins->mesh; 
  elliptic_t *solver = ins->hSolver; 
  
  if (ins->hOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    // ins->heatRhsBCKernel(mesh->Nelements,
    //                     mesh->o_ggeo,
    //                     mesh->o_sgeo,
    //                     mesh->o_Dmatrices,
    //                     mesh->o_Smatrices,
    //                     mesh->o_MM,
    //                     mesh->o_vmapM,
    //                     mesh->o_EToB,
    //                     mesh->o_sMT,
    //                     ins->lambdaHeat,
    //                     time,
    //                     mesh->o_x,
    //                     mesh->o_y,
    //                     mesh->o_z,
    //                     ins->o_VmapB,
    //                     ins->o_rhsT);
    
  } else if (ins->hOptions.compareArgs("DISCRETIZATION","IPDG")) {
    occaTimerTic(mesh->device,"heatRhsIpdg");    
    ins->heatRhsIpdgBCKernel(mesh->Nelements,
                            mesh->o_vmapM,
                            solver->tau,
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
                            ins->o_rhsT);
    occaTimerToc(mesh->device,"heatRhsIpdg");   
  }

// gather-scatter
  if(ins->hOptions.compareArgs("DISCRETIZATION","CONTINUOUS")){
    ogsGatherScatter(ins->o_rhsT, ogsDfloat, ogsAdd, mesh->ogs);
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_rhsT);
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_rkT);
  }

  occaTimerTic(mesh->device,"Heat Solve");
  ins->NiterT = ellipticSolve(solver, ins->lambdaHeat, ins->heatTOL, ins->o_rhsT, ins->o_rkT); 
  occaTimerToc(mesh->device,"Heat Solve"); 


  if (ins->hOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    ins->heatAddBCKernel(mesh->Nelements,
                        time,
                        mesh->o_sgeo,
                        mesh->o_x,
                        mesh->o_y,
                        mesh->o_z,
                        mesh->o_vmapM,
                        ins->o_VmapB,
                        ins->o_rkT);
  }
}
