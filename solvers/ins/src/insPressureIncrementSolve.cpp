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

#include "ins.hpp"

//  Solves -gamma*Laplacian*PI = rhs
//  P += PI
void ins_t::PressureIncrementSolve(occa::memory& o_P, occa::memory& o_RHS,
                                   const dfloat gamma, const dfloat T, const dfloat dt){

  // compute RHS = MM*RHS/gamma + BCdata
  pressureIncrementRhsKernel(mesh.Nelements,
                    mesh.o_vgeo,
                    mesh.o_sgeo,
                    mesh.o_ggeo,
                    mesh.o_S,
                    mesh.o_D,
                    mesh.o_LIFT,
                    mesh.o_MM,
                    mesh.o_sM,
                    mesh.o_vmapM,
                    mesh.o_EToB,
                    o_mapB,
                    pTau,
                    T,
                    dt,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    nu,
                    gamma,
                    o_RHS);

  // gather-scatter
  if(pDisc_c0) {
    mesh.ogs->GatherScatter(o_RHS, ogs_dfloat, ogs_add, ogs_sym);

    //PI should always be C0 so skip the gs+weighting

    if (pSolver->Nmasked) pSolver->maskKernel(pSolver->Nmasked, pSolver->o_maskIds, o_P);
    if (pSolver->Nmasked) pSolver->maskKernel(pSolver->Nmasked, pSolver->o_maskIds, o_RHS);
  }

  int maxIter = 5000;
  int verbose = 0;

  //  Solve - Laplacian*PI = RHS
  NiterP = pSolver->Solve(*pLinearSolver, o_PI, o_RHS, presTOL, maxIter, verbose);

  // P += PI and enter BCs if C0
  pressureIncrementBCKernel(mesh.Nelements,
                   mesh.o_sgeo,
                   mesh.o_vmapM,
                   o_mapB,
                   T,
                   mesh.o_x,
                   mesh.o_y,
                   mesh.o_z,
                   pDisc_c0,
                   nu,
                   o_PI,
                   o_P);
}
