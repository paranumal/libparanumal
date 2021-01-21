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

#include "ellipticPrecon.hpp"

// Jacobi preconditioner
JacobiPrecon::JacobiPrecon(elliptic_t& _elliptic):
  elliptic(_elliptic) {

  dfloat *diagA    = (dfloat*) calloc(elliptic.Ndofs, sizeof(dfloat));
  dfloat *invDiagA = (dfloat*) calloc(elliptic.Ndofs, sizeof(dfloat));
  elliptic.BuildOperatorDiagonal(diagA);
  for (dlong n=0;n<elliptic.Ndofs;n++)
    invDiagA[n] = 1.0/diagA[n];

  o_invDiagA = elliptic.platform.malloc(elliptic.Ndofs*sizeof(dfloat), invDiagA);

  free(diagA);
  free(invDiagA);
}

void JacobiPrecon::Operator(occa::memory& o_r, occa::memory& o_Mr) {

  // Mr = invDiag.*r
  elliptic.linAlg.amxpy(elliptic.Ndofs, 1.0, o_invDiagA, o_r, 0.0, o_Mr);

  // zero mean of RHS
  if(elliptic.allNeumann) elliptic.ZeroMean(o_Mr);
}