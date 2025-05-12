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

// 1. G*ML*S*GU = G*ML*U
// 2. U = S*U
void ins_t::MassSolve(deviceMemory<dfloat>& o_U){

  int Nfields = mesh.dim;
  
  dlong Ntotal = (mesh.Nelements+mesh.totalHaloPairs)*mesh.Np*Nfields;
  dlong Nglobal = (massSolver.Ndofs+massSolver.Nhalo);
  int maxIter = 1000;
  bool verbose = false;

  deviceMemory<dfloat> o_GrhsU = platform.reserve<dfloat>(Nglobal);
  deviceMemory<dfloat> o_GUH   = platform.reserve<dfloat>(Nglobal);
  deviceMemory<dfloat> o_rhsU  = platform.reserve<dfloat>(Ntotal);

  // compute RHS = MM*RHS/nu + BCdata
  // and split fields to separate arrays
  
  massSolver.massAxKernel(mesh.Nelements,
			  mesh.o_wJ,
			  mesh.o_MM,
			  o_U,
			  o_rhsU);
  
  // gather
  massSolver.ogsMasked.Gather(o_GrhsU, o_rhsU, Nfields, ogs::Add, ogs::Trans);

  // solve
  NiterU = massSolver.Solve(massLinearSolver, o_GUH, o_GrhsU, velTOL, maxIter, verbose);

  // scatter
  //  massSolver.ogsMasked.Scatter(o_U, o_GUH, Nfields1, ogs::NoTrans);

  massSolver.massScatterKernel(mesh.Nelements, massSolver.o_GlobalToLocal, o_GUH, o_U);
  
  o_GrhsU.free();
  o_GUH.free();
  o_rhsU.free();
  
}
