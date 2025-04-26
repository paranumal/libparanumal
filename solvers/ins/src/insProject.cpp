
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

void ins_t::Project(deviceMemory<dfloat>& o_U, int Nfilt){

  // weight and extract velocity

  dlong Nlocal = mesh.Nelements*mesh.Np;
  dlong Nhalo  = mesh.totalHaloPairs*mesh.Np;
  
  // scatter
  deviceMemory<dfloat> o_wUL = platform.reserve<dfloat>(Nlocal+Nhalo);
  deviceMemory<dfloat> o_wVL = platform.reserve<dfloat>(Nlocal+Nhalo);
  deviceMemory<dfloat> o_UG = platform.reserve<dfloat>(uSolver.Ndofs+uSolver.Nhalo);
  deviceMemory<dfloat> o_VG = platform.reserve<dfloat>(vSolver.Ndofs+vSolver.Nhalo);

  projectWeightKernel(mesh.Nelements, Nfilt, o_projectWeights, o_U, o_wUL, o_wVL);
  
  uSolver.ogsMasked.Gather(o_UG, o_wUL, 1, ogs::Add, ogs::Trans);
  if(Nfilt>1)
    vSolver.ogsMasked.Gather(o_VG, o_wVL, 1, ogs::Add, ogs::Trans);

  projectScatterKernel(mesh.Nelements, Nfilt, o_uGlobalToLocal, o_UG, o_vGlobalToLocal, o_VG, o_U);
  
}
