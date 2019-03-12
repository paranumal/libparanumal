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

#include "cds.h"

void cdsHelmholtzRhs(cds_t *cds, dfloat time, int stage, occa::memory o_rhsS){
  
  mesh_t *mesh = cds->mesh; 

  // rhsU^s = MM*(\sum^s b_i U^n-i - \sum^s-1 a_i N(U^n-i) + \sum^s-1 c_i GP^n-i)/nu dt
  cds->helmholtzRhsKernel(mesh->Nelements,
                          mesh->o_vgeo,
                          mesh->o_MM,
                          cds->idt,
                          cds->ialf,
                          cds->o_extbdfA,
                          cds->o_extbdfB,
                          cds->o_extbdfC,
                          cds->sOffset,
                          cds->o_S,
                          cds->o_NS,
                          o_rhsS);
}
