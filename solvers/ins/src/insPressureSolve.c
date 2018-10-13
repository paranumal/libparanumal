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

// solve Ap = rhsP
void insPressureSolve(ins_t *ins, dfloat time, int stage){

  mesh_t *mesh = ins->mesh;
  elliptic_t *solver = ins->pSolver;

  int quad3D = (ins->dim==3 && ins->elementType==QUADRILATERALS) ? 1 : 0;  

  if (ins->pOptions.compareArgs("DISCRETIZATION","CONTINUOUS") && !quad3D) {
    ins->pressureRhsBCKernel(mesh->Nelements,
                            mesh->o_ggeo,
                            mesh->o_sgeo,
                            mesh->o_Dmatrices,
                            mesh->o_Smatrices,
                            mesh->o_vmapM,
                            mesh->o_EToB,
                            mesh->o_sMT,
                            time,
                            ins->dt,
                            stage,
                            ins->ARKswitch,
                            ins->o_rkC,
                            ins->o_prkB,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            ins->o_PmapB,
                            ins->o_rhsP);
  } else if(ins->pOptions.compareArgs("DISCRETIZATION","IPDG") && !quad3D) {
    occaTimerTic(mesh->device,"PoissonRhsIpdg"); 
    ins->pressureRhsIpdgBCKernel(mesh->Nelements,
                                  mesh->o_vmapM,
                                  solver->tau,
                                  time,
                                  ins->dt,
                                  stage,
                                  ins->ARKswitch,
                                  ins->o_rkC,
                                  ins->o_prkB,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_Dmatrices,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  ins->o_rhsP);
    occaTimerToc(mesh->device,"PoissonRhsIpdg");
  }

  //keep current PI as the initial guess?

  // gather-scatter
  if(ins->pOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    ogsGatherScatter(ins->o_rhsP, ogsDfloat, ogsAdd, mesh->ogs);
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_rhsP);
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, ins->o_PI);
  }

  occaTimerTic(mesh->device,"Pr Solve");
  ins->NiterP = ellipticSolve(solver, 0.0, ins->presTOL, ins->o_rhsP, ins->o_PI); 
  occaTimerToc(mesh->device,"Pr Solve"); 

 if (ins->pOptions.compareArgs("DISCRETIZATION","CONTINUOUS") && !quad3D) {
    ins->pressureAddBCKernel(mesh->Nelements,
                            time,
                            ins->dt,
                            stage,
                            ins->ARKswitch,
                            ins->o_rkC,
                            ins->o_prkB,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            ins->o_PmapB,
                            ins->o_PI);
  }
}
