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

#include "ins.h"

// compute filter term and add to explicit storage data FU
void insAddVelocityRhs(ins_t *ins, dfloat time){

  mesh_t *mesh = ins->mesh;
  
  // Explicitly make FU zero to prevent multiple additions
  dfloat zero = 0.0; 
  ins->setScalarKernel(ins->Ntotal*ins->NVfields, zero, ins->o_FU); 
   
  if(ins->options.compareArgs("FILTER STABILIZATION", "RELAXATION"))
  ins->filterKernel(mesh->Nelements,
                    ins->o_filterMT,
                    ins->filterS, 
                    ins->fieldOffset,
                    ins->o_U,
                    ins->o_FU);
}
