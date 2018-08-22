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

// compute NU = N(U)
void mppfAdvection(mppf_t *mppf, dfloat time){

  mesh_t *mesh = mppf->mesh;
  
  // Compute Volume Contribution
  occaTimerTic(mesh->device,"AdvectionVolume");
  if(mppf->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
    mppf->advectionCubatureVolumeKernel(mesh->Nelements,
                                        mesh->o_vgeo,
                                        mesh->o_cubvgeo,
                                        mesh->o_cubDWmatrices,
                                        mesh->o_cubInterpT,
                                        mesh->o_cubProjectT,
                                        mppf->fieldOffset,
                                        mppf->o_U,
                                        mppf->o_cU,
                                        mppf->o_NU);
  } else {
    mppf->advectionVolumeKernel(mesh->Nelements,
                                mesh->o_vgeo,
                                mesh->o_Dmatrices,
                                mppf->fieldOffset,
                                mppf->o_U,
                                mppf->o_NU);
  }
  occaTimerToc(mesh->device,"AdvectionVolume");

 
  occaTimerTic(mesh->device,"AdvectionSurface");
  if(mppf->options.compareArgs("ADVECTION TYPE", "CUBATURE")){
    mppf->advectionCubatureSurfaceKernel(mesh->Nelements,
                                          mesh->o_vgeo,
                                          mesh->o_sgeo,
                                          mesh->o_cubsgeo,
                                          mesh->o_intInterpT,
                                          mesh->o_intLIFTT,
                                          mesh->o_cubInterpT,
                                          mesh->o_cubProjectT,
                                          mesh->o_vmapM,
                                          mesh->o_vmapP,
                                          mesh->o_EToB,
                                          time,
                                          mesh->o_intx,
                                          mesh->o_inty,
                                          mesh->o_intz,
                                          mppf->fieldOffset,
                                          mppf->o_U,
                                          mppf->o_NU);
  } else {
    mppf->advectionSurfaceKernel(mesh->Nelements,
                                  mesh->o_sgeo,
                                  mesh->o_LIFTT,
                                  mesh->o_vmapM,
                                  mesh->o_vmapP,
                                  mesh->o_EToB,
                                  time,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_z,
                                  mppf->fieldOffset,
                                  mppf->o_U,
                                  mppf->o_NU);
  }
  occaTimerToc(mesh->device,"AdvectionSurface");


// Give exact NU
#if 0
  for(int e=0; e<mesh->Nelements;e++){
    for(int n=0; n<mesh->Np; n++){
      const int id = e*mesh->Np + n;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];

      dfloat nux = -(M_PI*sin(2*M_PI*x)*(cos(2.0*time)/2.0 - 0.5))/2.0;
      dfloat nuy = -(M_PI*sin(2*M_PI*y)*(cos(2.0*time)/2.0 - 0.5))/2.0;

      mppf->rkU[id + 0*mppf->fieldOffset] = nux;
      mppf->rkU[id + 1*mppf->fieldOffset] = nuy;
    }
  }
  
  mppf->o_GSave.copyFrom(mppf->rkU, mppf->NVfields*mppf->Ntotal*sizeof(dfloat));
  o_NU.copyFrom(mppf->o_GSave, mppf->NVfields*mppf->Ntotal*sizeof(dfloat), 0, 0); 

#endif






}
