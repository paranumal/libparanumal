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

void adaptiveOperator(adaptive_t *adaptive,
		      level_t *level,
		      dfloat lambda,
		      occa::memory &o_q,
		      occa::memory &o_Aq){

  setupAide &options = adaptive->options;

  int enableGatherScatters = 1;
  int enableReductions = 1;
  int continuous = options.compareArgs("DISCRETIZATION", "CONTINUOUS");
  int serial = options.compareArgs("THREAD MODEL", "Serial");
  int ipdg = options.compareArgs("DISCRETIZATION", "IPDG");
  
  if(continuous){

    ogs_t *ogs = level->ogs;
    
    dfloat_t lambda = 1.0;
    level->compute_partial_Ax(level->Klocal, // locally owned elements
			      level->o_IToE,
			      level->o_ggeo,
			      level->o_D,
			      lambda,
			      o_q,
			      o_Aq);

    adaptiveGatherScatter(adaptive, level, o_Aq);

    // boost null space option
#if USE_NULL_BOOST==1
    if(level->allNeumann) {

      dfloat alpha = 0., alphaG = 0.;
      dlong Nblock = level->Nblock;
      dfloat *tmp = level->tmp;
      occa::memory &o_tmp = level->o_tmp;
      
      adaptive->innerProductKernel(level->Nelements*level->Np, level->o_invDegree, o_q, o_tmp);
      o_tmp.copyTo(tmp);
      
      for(dlong n=0;n<Nblock;++n)
        alpha += tmp[n];

      MPI_Allreduce(&alpha, &alphaG, 1, MPI_DFLOAT, MPI_SUM, level->comm);
      alphaG *= level->allNeumannPenalty*level->allNeumannScale*level->allNeumannScale;

      adaptive->addScalarKernel(level->Nelements*level->Np, alphaG, o_Aq);
    }
#endif

#if 0
    //post-mask
    if (level->Nmasked) 
      adaptive->maskKernel(level->Nmasked, level->o_maskIds, o_Aq);
#endif

  } else if(ipdg){
    printf("WARNING: bailing since ipdg is not yet supported\n");
    MPI_Finalize();
    exit(0);
  } 

}

