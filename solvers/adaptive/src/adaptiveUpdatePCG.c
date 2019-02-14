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


dfloat adaptiveUpdatePCG(adaptive_t *adaptive,
			 occa::memory &o_p,
			 occa::memory &o_Ap,
			 const dfloat alpha,
			 occa::memory &o_x,
			 occa::memory &o_r){

  setupAide &options = adaptive->options;
  
  int fixedIterationCountFlag = 0;
  int enableGatherScatters = 1;
  int enableReductions = 1;
  int flexible = options.compareArgs("KRYLOV SOLVER", "FLEXIBLE");
  int verbose = options.compareArgs("VERBOSE", "TRUE");
  int serial = options.compareArgs("THREAD MODEL", "Serial");
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int ipdg = options.compareArgs("DISCRETIZATION", "IPDG");
  

  int Ntotal = level->Klocal*level->Np;
  
  dfloat rdotr1 = 0;
  
  if(continuous){
    
    // x <= x + alpha*p
    // r <= r - alpha*A*p
    // dot(r,r)
    adaptive->updatePCGKernel(Ntotal,
			      level->NblocksUpdatePCG,
			      level->o_invDegree, o_p, o_Ap, alpha, o_x, o_r, level->o_tmpNormr);

    level->o_tmpNormr.copyTo(level->tmpNormr);
    
    rdotr1 = 0;
    for(int n=0;n<level->NblocksUpdatePCG;++n){
      rdotr1 += level->tmpNormr[n];
    }
    
    dfloat globalrdotr1 = 0;
    MPI_Allreduce(&rdotr1, &globalrdotr1, 1, MPI_DFLOAT, MPI_SUM, mesh->comm);

    rdotr1 = globalrdotr1;
    
  }

  return rdotr1;
}

