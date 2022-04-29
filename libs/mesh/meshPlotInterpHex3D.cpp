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

#include "mesh.hpp"

namespace libp {

//interpolate field to plotting nodes
void mesh_t::PlotInterpHex3D(const memory<dfloat> q, memory<dfloat> Iq, memory<dfloat> scratch){

  if (scratch.length()< static_cast<size_t>(plotNq*Nq*Nq + plotNq*plotNq*Nq)) {
    //if not provided with enough scratch space, alloc our own
    scratch.malloc(plotNq*Nq*Nq + plotNq*plotNq*Nq);
  }

  memory<dfloat> IQ  = scratch;
  memory<dfloat> IIQ = scratch + plotNq*Nq*Nq;

  //interpolate in r
  for(int k=0;k<Nq;++k){
    for(int j=0;j<Nq;++j){
      for(int i=0;i<plotNq;++i){
        dfloat qn = 0;

        for(int m=0;m<Nq;++m){
          const int qid = m + j*Nq + k*Nq*Nq;
          qn += plotInterp[i*Nq+m]*q[qid];
        }

        const int id = i + j*plotNq + k*Nq*plotNq;
        IQ[id] = qn;
      }
    }
  }

  //interpolate in s
  for(int k=0;k<Nq;++k){
    for(int j=0;j<plotNq;++j){
      for(int i=0;i<plotNq;++i){
        dfloat qn = 0;

        for(int m=0;m<Nq;++m){
          const int qid = i + m*plotNq + k*Nq*plotNq;
          qn += plotInterp[j*Nq+m]*IQ[qid];
        }

        const int id = i + j*plotNq + k*plotNq*plotNq;
        IIQ[id] = qn;
      }
    }
  }

  //interpolate in k and write
  for(int k=0;k<plotNq;++k){
    for(int j=0;j<plotNq;++j){
      for(int i=0;i<plotNq;++i){
        dfloat qn = 0;

        for(int m=0;m<Nq;++m){
          const int qid = i + j*plotNq + m*plotNq*plotNq;
          qn += plotInterp[k*Nq+m]*IIQ[qid];
        }

        const int id = i + j*plotNq + k*plotNq*plotNq;
        Iq[id] = qn;
      }
    }
  }
}

} //namespace libp
