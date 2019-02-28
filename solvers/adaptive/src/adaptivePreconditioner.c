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
                            occa::memory &o_r, occa::memory &o_tmp, occa::memory &o_z){

  level_t *level0 = adaptive->lvl;
  setupAide options = adaptive->options;

  dlong Ntotal = level0->Np*level0->Klocal;

  if(options.compareArgs("PRECONDITIONER", "JACOBI")){

    // Jacobi preconditioner     (builds in W)
    adaptive->dotMultiplyKernel(Ntotal, level0->o_invDiagA,  o_r, o_z);

    adaptive->dotMultiplyKernel(Ntotal, level0->o_invDegree, o_z, o_z);

    ogsGatherScatter(o_z, ogsDfloat, ogsAdd, level0->ogs);
    //    adaptiveGatherScatter(level0, o_z);
    
#if 1    
#if USE_GASPAR==1
    level0->scatter_noncon(level0->Klocal, level0->o_EToC,
			   level0->o_Ib, level0->o_It, o_z);
#endif 
#endif
    
    //    adaptiveRankOneProjection(level0, level0->o_filtU, level0->o_filtV, o_z, o_z);
    
  }
  else if (options.compareArgs("PRECONDITIONER", "MULTIGRID")) {

    for(int L=0;L<adaptive->Nlevels;++L){

      // smooth at level L
      
      // restrict from level L to coarser grid (L+1)

    }

    // smooth at level L = adaptive->Nlevels-1
    
    for(int L=adaptive->Nlevels;L>=0;--L){

      // prolongate from level L to finer grid (L-1)
      
      // smooth at level (L-1)
    }
    
  }
  else{ // turn off preconditioner

#if USE_GASPAR==1
    adaptive->dotMultiplyKernel(Ntotal, level0->o_invDegree, o_r, o_z);
    
    ogsGatherScatter(o_z, ogsDfloat, ogsAdd, level0->ogs);
    
    level0->scatter_noncon(level0->Klocal, level0->o_EToC,
			   level0->o_Ib, level0->o_It, o_z);
#else
    adaptive->dotMultiplyKernel(Ntotal, level0->o_invDegree, o_r, o_z);
    
    ogsGatherScatter(o_z, ogsDfloat, ogsAdd, level0->ogs);

#endif

    o_z.copyFrom(o_r);
    
  }
  
#if USE_NULL_PROJECTION==1
  if(adaptive->allNeumann) // zero mean of RHS
    adaptiveZeroMean(adaptive, level0, o_z);
#endif

}


