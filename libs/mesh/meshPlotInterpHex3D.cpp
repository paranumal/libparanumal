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

#include "mesh.hpp"
#include "mesh/mesh3D.hpp"

//interpolate field to plotting nodes
void meshHex3D::PlotInterp(const dfloat* q, dfloat* Iq, dfloat* scratch){

  dfloat *IQ, *IIQ;

  bool alloc_scratch=false;
  if (scratch==nullptr) {
    //if not provided with a scratch space, alloc our own
    alloc_scratch=true;
    IQ  = (dfloat *) malloc(plotNq*Nq*Nq*sizeof(dfloat));
    IIQ = (dfloat *) malloc(plotNq*plotNq*Nq*sizeof(dfloat));
  } else {
    IQ  = scratch;
    IIQ = scratch + plotNq*Nq*Nq;
  }

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

  //clean up
  if (alloc_scratch) {
    free(IQ); free(IIQ);
  }
}
