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

//  Solves gamma*U - nu*Laplacian*U = rhs
void ins_t::VelocitySolve(deviceMemory<dfloat>& o_U, deviceMemory<dfloat>& o_RHS,
                          const dfloat gamma, const dfloat T) {

  // compute RHS = MM*RHS/nu + BCdata
  // and split fields to separate arrays
  velocityRhsKernel(mesh.Nelements,
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
                      vTau,
                      T,
                      mesh.o_x,
                      mesh.o_y,
                      mesh.o_z,
                      gamma/nu,
                      nu,
                      o_U,
                      o_RHS,
                      o_UH,
                      o_VH,
                      o_WH,
                      o_rhsU,
                      o_rhsV,
                      o_rhsW);

  int maxIter = 5000;
  int verbose = 0;

  uSolver.lambda = gamma/nu;
  vSolver.lambda = gamma/nu;
  wSolver.lambda = gamma/nu;

  //  Solve lambda*U - Laplacian*U = rhs
  if (vDisc_c0){
    // gather, solve, scatter
    uSolver.ogsMasked.Gather(o_GrhsU, o_rhsU, 1, ogs::Add, ogs::Trans);
    NiterU = uSolver.Solve(uLinearSolver, o_GUH, o_GrhsU, velTOL, maxIter, verbose);
    uSolver.ogsMasked.Scatter(o_UH, o_GUH, 1, ogs::NoTrans);

    vSolver.ogsMasked.Gather(o_GrhsV, o_rhsV, 1, ogs::Add, ogs::Trans);
    NiterV = vSolver.Solve(vLinearSolver, o_GVH, o_GrhsV, velTOL, maxIter, verbose);
    vSolver.ogsMasked.Scatter(o_VH, o_GVH, 1, ogs::NoTrans);

    if (mesh.dim==3) {
      wSolver.ogsMasked.Gather(o_GrhsW, o_rhsW, 1, ogs::Add, ogs::Trans);
      NiterW = wSolver.Solve(wLinearSolver, o_GWH, o_GrhsW, velTOL, maxIter, verbose);
      wSolver.ogsMasked.Scatter(o_WH, o_GWH, 1, ogs::NoTrans);
    }

  } else {
    NiterU = uSolver.Solve(uLinearSolver, o_UH, o_rhsU, velTOL, maxIter, verbose);
    NiterV = vSolver.Solve(vLinearSolver, o_VH, o_rhsV, velTOL, maxIter, verbose);
    if (mesh.dim==3)
      NiterW = wSolver.Solve(wLinearSolver, o_WH, o_rhsW, velTOL, maxIter, verbose);
  }

  // merge arrays back, and enter BCs if C0
  velocityBCKernel(mesh.Nelements,
                  mesh.o_sgeo,
                  mesh.o_vmapM,
                  mesh.o_mapB,
                  T,
                  mesh.o_x,
                  mesh.o_y,
                  mesh.o_z,
                  nu,
                  vDisc_c0,
                  o_UH,
                  o_VH,
                  o_WH,
                  o_U);
}
