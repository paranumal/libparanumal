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

void adaptivePreconditioner(adaptive_t *adaptive, dfloat lambda,
                            occa::memory &o_r, occa::memory &o_z){

  level_t *level0 = adaptive->lvl;
  setupAide options = adaptive->options;

  dlong Ntotal = level0->Np*level0->Klocal;

  
  if(options.compareArgs("PRECONDITIONER", "JACOBI")){
    
    // Jacobi preconditioner
    
    adaptive->dotMultiplyKernel(Ntotal, o_r, level0->o_invDiagA, o_z);
  }
  else if (options.compareArgs("PRECONDITIONER", "MULTIGRID")) {

  }
  else{ // turn off preconditioner

    adaptive->dotMultiplyKernel(Ntotal, level0->o_invDegree, o_r, o_z);
    
    ogsGatherScatter(o_z, ogsDfloat, ogsAdd, level0->ogs);
    
    level0->scatter_noncon(level0->Klocal, level0->o_EToC, level0->o_Ib, level0->o_It, o_z);

  }

#if USE_NULL_PROJECTION==1
  if(adaptive->allNeumann) // zero mean of RHS
    adaptiveZeroMean(adaptive, level0, o_z);
#endif
}


