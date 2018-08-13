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

void mppfPressureSolve(mppf_t *mppf, dfloat time, occa::memory o_rkP){

  mesh_t *mesh = mppf->mesh;
  elliptic_t *solver = mppf->pSolver;

  if (mppf->pOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    mppf->pressureRhsBCKernel(mesh->Nelements,
                            mesh->o_ggeo,
                            mesh->o_sgeo,
                            mesh->o_Dmatrices,
                            mesh->o_Smatrices,
                            mesh->o_vmapM,
			                      mesh->o_EToB,
                            mesh->o_sMT,
                            time,
                            mppf->dt,
                            mppf->Nstages,
                            mppf->ARKswitch,
                            mppf->o_rkC,
                            mppf->o_prkB,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mppf->o_PmapB,
                            mppf->o_rhsP);
  } else if (mppf->pOptions.compareArgs("DISCRETIZATION","IPDG")) {
    occaTimerTic(mesh->device,"PoissonRhsIpdg"); 
    mppf->pressureRhsIpdgBCKernel(mesh->Nelements,
                                  mesh->o_vmapM,
                                  solver->tau,
                                  time,
                                  mppf->dt,
                                  mppf->Nstages,
                                  mppf->ARKswitch,
                                  mppf->o_rkC,
                                  mppf->o_prkB,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_Dmatrices,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  mppf->o_rhsP);
    occaTimerToc(mesh->device,"PoissonRhsIpdg");
  }

  // gather-scatter
  if(mppf->pOptions.compareArgs("DISCRETIZATION","CONTINUOUS")){
    ellipticParallelGatherScatter(mesh, mesh->ogs, mppf->o_rhsP, dfloatString, "add");  
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, mppf->o_rhsP);
    if (solver->Nmasked) mesh->maskKernel(solver->Nmasked, solver->o_maskIds, mppf->o_rkP);
  }

  occaTimerTic(mesh->device,"Pr Solve");
  mppf->NiterP = ellipticSolve(solver, 0.0, mppf->presTOL, mppf->o_rhsP, mppf->o_rkP); 
  occaTimerToc(mesh->device,"Pr Solve"); 

  if (mppf->pOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
    mppf->pressureAddBCKernel(mesh->Nelements,
                            time,
                            mppf->dt,
                            mppf->Nstages,
                            mppf->ARKswitch,
                            mppf->o_rkC,
                            mppf->o_prkB,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_z,
                            mesh->o_vmapM,
                            mppf->o_PmapB,
                            mppf->o_rkP);
  }
}
