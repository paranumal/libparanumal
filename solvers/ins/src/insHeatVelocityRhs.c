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

void insHeatVelocityRhs(ins_t *ins, dfloat time, int stage, occa::memory o_rhsU, occa::memory o_rhsV, 
                                                            occa::memory o_rhsW, occa::memory o_rhsT){
  
  mesh_t *mesh = ins->mesh; 

  if (ins->options.compareArgs("TIME INTEGRATOR", "ARK")) { 
    printf("ARK is not implemented for heat solver\n");
    // // rhsU^s = MM*(U^n - \sum^s-1 ea_si N(U^i) + \sum^s-1 ia_si LU^i - \sum^s-1 pa_si GP^i)/ia_ss nu dt
    // ins->heatRhsKernel(mesh->Nelements,
    //                        stage,
    //                        mesh->o_vgeo,
    //                        mesh->o_MM,
    //                        ins->idt,
    //                        ins->ialpha,
    //                        ins->o_erkA,
    //                        ins->o_irkA,
    //                        ins->o_prkA,
    //                        ins->o_prkB,
    //                        ins->fieldOffset,
    //                        ins->o_U,
    //                        ins->o_T,
    //                        ins->o_NT,
    //                        ins->o_LU,
    //                        ins->o_GP,
    //                        o_rhsU,
    //                        o_rhsV,
    //                        o_rhsW);
  } else if (ins->options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    // rhsU^s = MM*(\sum^s b_i U^n-i - \sum^s-1 a_i N(U^n-i) + \sum^s-1 c_i GP^n-i)/nu dt
    ins->heatVelocityRhsKernel(mesh->Nelements,
                               mesh->o_vgeo,
                               mesh->o_MM,
                               ins->idt,
                               ins->inu,
                               ins->ialpha,
                               ins->o_extbdfA,
                               ins->o_extbdfB,
                               ins->o_extbdfC,
                               ins->fieldOffset,
                               ins->o_U,
                               ins->o_T,
                               ins->o_NU,
                               ins->o_NT,
                               ins->o_GP,
                               o_rhsU,
                               o_rhsV,
                               o_rhsW,
                               o_rhsT);
  }
}
