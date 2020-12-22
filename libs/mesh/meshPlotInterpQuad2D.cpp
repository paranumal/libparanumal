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
#include "mesh/mesh2D.hpp"

//interpolate field to plotting nodes
void meshQuad2D::PlotInterp(const dfloat* q, dfloat* Iq, dfloat* scratch){

  dfloat *IQ;

  bool alloc_scratch=false;
  if (scratch==nullptr) {
    //if not provided with a scratch space, alloc our own
    alloc_scratch=true;
    IQ  = (dfloat *) malloc(plotNq*Nq*sizeof(dfloat));
  } else {
    IQ  = scratch;
  }

  //interpolate in r
  for(int j=0;j<Nq;++j){
    for(int i=0;i<plotNq;++i){
      dfloat qn = 0;

      for(int m=0;m<Nq;++m){
        const int qid = m + j*Nq;
        qn += plotInterp[i*Nq+m]*q[qid];
      }

      const int id = i + j*plotNq;
      IQ[id] = qn;
    }
  }

  //interpolate in s
  for(int j=0;j<plotNq;++j){
    for(int i=0;i<plotNq;++i){
      dfloat qn = 0;

      for(int m=0;m<Nq;++m){
        const int qid = i + m*plotNq;
        qn += plotInterp[j*Nq+m]*IQ[qid];
      }

      const int id = i + j*plotNq;
      Iq[id] = qn;
    }
  }

  //clean up
  if (alloc_scratch) {
    free(IQ);
  }
}
