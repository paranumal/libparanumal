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

#include "mppf.h"

void mppfPressureRhs(mppf_t *mppf, dfloat time, occa::memory o_rkU){

  mesh_t *mesh = mppf->mesh;


  // rhsP = Div Uhat
  mppfDivergence(mppf, time, o_rkU, mppf->o_rhsP);
  
  // // rhsP = -MM*rho0* Div Uhat * invdt
  occaTimerTic(mesh->device,"PoissonRhsForcing");
  mppf->pressureRhsKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_MM,
                              mppf->idt,  
                              mppf->o_rhsP);
  occaTimerToc(mesh->device,"PoissonRhsForcing");

  
  

#if 0
  //add penalty from jumps in previous pressure
  mppf->poissonPenaltyKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_vgeo,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                mppf->tau,
                                mesh->o_x,
                                mesh->o_y,
                                t,
                                mppf->dt,
                                mppf->c0,
                                mppf->c1,
                                mppf->c2,
                                mppf->index,
                                (mesh->Nelements+mesh->totalHaloPairs),
                                mppf->o_P,
                                mppf->o_rhsP);
  #endif
}
