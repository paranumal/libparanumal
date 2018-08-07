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

// complete a time step using LSERK4
void mppfAdvectionUpdate(mppf_t *mppf, dfloat time, occa::memory o_NU, occa::memory o_GU, occa::memory o_GP, occa::memory o_rkU){

  mesh_t *mesh = mppf->mesh;
  
  mppf->advectionUpdateKernel(mesh->Nelements,
                           mesh->o_vgeo,
                           mesh->o_Dmatrices,
                           mesh->o_x,
                           mesh->o_y,
                           mesh->o_z,
                           mppf->dt,
                           mppf->time,
                           mppf->o_extbdfA,
                           mppf->o_extbdfB,
                           mppf->fieldOffset,
                           mppf->o_U,
                           mppf->o_Rho,
                           mppf->o_Mu,
                           mppf->o_GMu,
                           o_NU,
                           o_GU,
                           o_GP,
                           mppf->o_Psi,
                           mppf->o_GPhi,
                           o_rkU);
}
