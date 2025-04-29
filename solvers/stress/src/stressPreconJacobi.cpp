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

  memory<dfloat> diagA   (stress.Ndofs);
  memory<pfloat> invDiagA(stress.Ndofs);
  stress.BuildOperatorDiagonal(diagA);

  for (dlong n=0;n<stress.Ndofs;n++){
    invDiagA[n] = 1.0/diagA[n];
  }

  o_invDiagA = stress.platform.malloc<pfloat>(invDiagA);

  std::cout << "Ndofs: " << stress.Ndofs << std::endl;
  
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
