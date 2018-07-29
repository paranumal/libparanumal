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

#include "mesh2D.h"

// function to compute collocation differentiation
// contributions to nodal DG rhs for acoustics
void acousticsPml2D(mesh2D *mesh){

  // for all elements
  for(int m=0;m<mesh->pmlNelements;++m){
    int e = mesh->pmlElementList[m];

    // prefetch geometric factors (constant on triangle)
    dfloat sigmax = mesh->pmlSigmaX[m];
    dfloat sigmay = mesh->pmlSigmaY[m];

    // for all nodes in this element
    for(int n=0;n<mesh->Np;++n){

      // load state and rhs values at node
      int   base = mesh->Nfields*(n+e*mesh->Np);
      dfloat u = mesh->q[base+0];
      dfloat v = mesh->q[base+1];
      dfloat p = mesh->q[base+2];

      // load pml state at node
      int   pmlBase = mesh->pmlNfields*(n+m*mesh->Np);
      dfloat utilde = mesh->pmlq[pmlBase+0];
      dfloat vtilde = mesh->pmlq[pmlBase+1];
      dfloat ptilde = mesh->pmlq[pmlBase+2];

      // update for u,v,p
      dfloat rhsu = -(sigmax-sigmay)*u - sigmay*(sigmay-sigmax)*utilde; // uhat
      dfloat rhsv = -(sigmay-sigmax)*v - sigmax*(sigmax-sigmay)*vtilde; // vhat
      dfloat rhsp = -(sigmax+sigmay)*p - sigmax*sigmay*ptilde; // p

      // update for u~,v~, p~
      dfloat rhsutilde = u-sigmay*utilde; // du~/dt = -sigmay*u~  + uhat
      dfloat rhsvtilde = v-sigmax*vtilde; // dv~/dt = -sigmax*v~  + vhat
      dfloat rhsptilde = p;                  // dp~/dt = p

      // store acoustics rhs contributions from PML
      mesh->rhsq[base+0] += rhsu;
      mesh->rhsq[base+1] += rhsv;
      mesh->rhsq[base+2] += rhsp;

      // store PML rhs
      mesh->pmlrhsq[pmlBase+0] = rhsutilde;
      mesh->pmlrhsq[pmlBase+1] = rhsvtilde;
      mesh->pmlrhsq[pmlBase+2] = rhsptilde;

    }
  }
}
