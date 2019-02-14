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

#include "adaptive.h"
#include "ogsInterface.h"

void adaptiveZeroMean(adaptive_t *adaptive, level_t *level, occa::memory &o_q){

  dfloat qmeanLocal;
  dfloat qmeanGlobal;
  
  dlong Nblock = level->Nblock;
  dfloat *tmp = level->tmp;

  occa::memory &o_tmp = level->o_tmp;

  // this is a C0 thing [ assume GS previously applied to o_q ]
#define USE_WEIGHTED 0
  
#if USE_WEIGHTED==1
  adaptive->innerProductKernel(level->Klocal*level->Np, level->o_invDegree, o_q, o_tmp);
#else
  adaptive->sumKernel(level->Klocal*level->Np, o_q, o_tmp);
#endif
  
  o_tmp.copyTo(tmp);

  // finish reduction
  qmeanLocal = 0;
  for(dlong n=0;n<Nblock;++n)
    qmeanLocal += tmp[n];

  // globalize reduction
  MPI_Allreduce(&qmeanLocal, &qmeanGlobal, 1, MPI_DFLOAT, MPI_SUM, adaptive->comm);

  // normalize
#if USE_WEIGHTED==1
  qmeanGlobal *= level->nullProjectWeightGlobal;
#else
  qmeanGlobal /= ((dfloat) level->Kglobal*(dfloat)level->Np);
#endif
  
  // q[n] = q[n] - qmeanGlobal
  adaptive->addScalarKernel(level->Klocal*level->Np, -qmeanGlobal, o_q);
}
