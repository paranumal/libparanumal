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
void mppfExplicitDiffusive(mppf_t *mppf, dfloat time, occa::memory o_U, occa::memory o_DU, occa::memory o_SU){

  mesh_t *mesh = mppf->mesh;
  // Halo exchange is already done in Advection kernel, just extrapolate elements with halo 
  const dlong NtotalElements = mesh->Nelements+mesh->totalHaloPairs;  
  
  occaTimerTic(mesh->device,"velocityExtrapolate");
  //compute advective velocity fields at time t
  mppf->velocityExtrapolateKernel(NtotalElements,
                                   mppf->Nstages,
                                   mppf->fieldOffset,
                                   mppf->o_extbdfA,
                                   o_U,
                                   mppf->o_Ue);

  occaTimerToc(mesh->device,"velocityExtrapolate");

// COMPLETE HALO EXCHANGE
  if (mppf->vOptions.compareArgs("DISCRETIZATION","CONTINUOUS")) {
  occaTimerTic(mesh->device,"velocityGradientVolume");
  // Compute Volume Contribution
  mppf->explicitDiffusiveKernel(mesh->Nelements,
                                      mesh->o_vgeo,
                                      mesh->o_Dmatrices,
                                      mppf->fieldOffset,
                                      mppf->o_Rho,
                                      mppf->o_Mu,
                                      mppf->o_Phi,
                                      mppf->o_Ue,
                                      o_DU,
                                      o_SU);
  occaTimerToc(mesh->device,"velocityGradientVolume");
}
else if (mppf->vOptions.compareArgs("DISCRETIZATION","IPDG")) {
  occaTimerTic(mesh->device,"velocityGradientVolume");
 
    occaTimerTic(mesh->device,"explicitDiffusiveVolume");
    // Compute Volume Contribution
    mppf->explicitDiffusiveVolumeKernel(mesh->Nelements,
                                      mesh->o_vgeo,
                                      mesh->o_Dmatrices,
                                      mppf->fieldOffset,
                                      mppf->o_Ue,
                                      mppf->o_GU);
    occaTimerToc(mesh->device,"explicitDiffusiveVolume");



    occaTimerTic(mesh->device,"explicitDiffusiveSurface");
    // Compute Surface Conribution
    mppf->explicitDiffusiveSurfaceKernel(mesh->Nelements,
                                         mesh->o_vgeo,
                                         mesh->o_sgeo,
                                         mesh->o_Dmatrices,
                                         mesh->o_LIFTT,
                                         mesh->o_vmapM,
                                         mesh->o_vmapP,
                                         mesh->o_EToB,
                                         mesh->o_x,
                                         mesh->o_y,
                                         mesh->o_z,
                                         mppf->o_Rho,
                                         mppf->o_Mu,
                                         mppf->o_GPhi,
                                         time,
                                         mppf->fieldOffset,
                                         mppf->o_Ue,
                                         mppf->o_GU,
                                         o_DU,
                                         o_SU);
    occaTimerToc(mesh->device,"explicitDiffusiveSurface");
  }
}
