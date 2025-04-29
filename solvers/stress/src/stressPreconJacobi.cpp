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

#include "stressPrecon.hpp"

// Jacobi preconditioner
JacobiPrecon::JacobiPrecon(stress_t& _stress):
  stress(_stress) {

#if 0
  memory<dfloat> diagA   (stress.Ndofs);
  memory<pfloat> invDiagA(stress.Ndofs);
  stress.BuildOperatorDiagonal(diagA);

  for (dlong n=0;n<stress.Ndofs;n++){
    invDiagA[n] = 1.0/diagA[n];
  }

  o_invDiagA = stress.platform.malloc<pfloat>(invDiagA);
#else
  
  deviceMemory<dfloat> o_diagAL =
    stress.platform.malloc<dfloat>(stress.mesh.Nelements*stress.Nfields*stress.mesh.Np);
  deviceMemory<dfloat> o_diagA =
    stress.platform.malloc<dfloat>(stress.Ndofs);

  o_invDiagA = stress.platform.malloc<pfloat>(stress.Ndofs);

  dfloat neumannBoost = (stress.allNeumann) ? stress.allNeumannPenalty*
    stress.allNeumannScale*stress.allNeumannScale: 0;

  stress.buildOperatorDiagonalKernel(stress.mesh.Nelements,
				     stress.o_nut,
				     stress.mesh.o_mapB,
				     neumannBoost,
				     stress.mesh.o_wJ,
				     stress.mesh.o_vgeo,
				     stress.mesh.o_D,
				     stress.mesh.o_S,
				     stress.mesh.o_MM,
				     stress.lambda,
				     o_diagAL);
  
  stress.ogsMasked.Gather(o_diagA, o_diagAL, stress.Nfields, ogs::Add, ogs::Trans);

  stress.reciprocalKernel(stress.Ndofs, o_diagA, o_invDiagA);

#endif
}

void JacobiPrecon::Operator(deviceMemory<pfloat>& o_r, deviceMemory<pfloat>& o_Mr) {

  linAlg_t& linAlg = stress.platform.linAlg();

  // Mr = invDiag.*r
  linAlg.amxpy(stress.Ndofs, (pfloat)1.0, o_invDiagA, o_r, (pfloat)0.0, o_Mr);

  // zero mean of RHS
  if(stress.allNeumann){
    stress.ZeroMean(o_Mr);
  }
}
