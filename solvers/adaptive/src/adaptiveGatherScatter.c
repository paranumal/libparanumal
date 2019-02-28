/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Ali Karakus, Lucas Wilcox

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

void adaptiveGatherScatter(level_t *level, occa::memory &o_x){

  //  adaptiveRankOneProjection(level, level->o_filtU, level->o_filtV, o_x, o_x);
  
  // gather over noncon faces to coarse side dofs
  level->gather_noncon(level->Klocal, level->o_EToC, level->o_Ib, level->o_It, o_x);
  
#ifdef PLUMADG_GATHER_SCATTER
  // FIXME for MPI
  level->gather_scatter(level->Ncontinuous, level->o_CToD_starts, level->o_CToD_indices, o_x);
#else
  ogsGatherScatter(o_x, ogsDfloat, ogsAdd, level->ogs);
#endif

  // scatter from coarse to fine noncon
  level->scatter_noncon(level->Klocal, level->o_EToC, level->o_Ib, level->o_It, o_x);

  //  adaptiveRankOneProjection(level, level->o_filtV, level->o_filtU, o_x, o_x);
  
}
