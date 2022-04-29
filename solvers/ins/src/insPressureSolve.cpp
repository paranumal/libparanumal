/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

//  Solves -gamma*Laplacian*P = rhs
void ins_t::PressureSolve(deviceMemory<dfloat>& o_P, deviceMemory<dfloat>& o_RHS,
                          const dfloat gamma, const dfloat T){

  // compute RHS = MM*RHS/gamma + BCdata
  pressureRhsKernel(mesh.Nelements,
                    mesh.o_wJ,
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
                    mesh.o_mapB,
                    pTau,
                    T,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    nu,
                    gamma,
                    o_RHS);

  //  Solve - Laplacian*P = RHS
  int maxIter = 5000;
  int verbose = 0;

  if(pDisc_c0) {
    // gather, solve, scatter
    pSolver.ogsMasked.Gather(o_GrhsP, o_RHS, 1, ogs::Add, ogs::Trans);
    NiterP = pSolver.Solve(pLinearSolver, o_GP, o_GrhsP, presTOL, maxIter, verbose);
    pSolver.ogsMasked.Scatter(o_P, o_GP, 1, ogs::NoTrans);

    // enter BCs if C0
    pressureBCKernel(mesh.Nelements,
                     mesh.o_sgeo,
                     mesh.o_vmapM,
                     mesh.o_mapB,
                     T,
                     mesh.o_x,
                     mesh.o_y,
                     mesh.o_z,
                     nu,
                     o_P);
  } else {
    NiterP = pSolver.Solve(pLinearSolver, o_P, o_RHS, presTOL, maxIter, verbose);
  }
}
