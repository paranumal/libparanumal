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

#include "bns.hpp"

//evaluate ODE rhs = f(q,t)
void bns_t::rhsf_pml(occa::memory& o_Q, occa::memory& o_pmlQ,
                     occa::memory& o_RHS, occa::memory& o_pmlRHS, const dfloat T){

  // extract q trace halo and start exchange
  traceHalo->ExchangeStart(o_Q, 1, ogs_dfloat);

  // compute volume contribution to bns RHS
  rhsVolume(mesh.NnonPmlElements, mesh.o_nonPmlElements, o_Q, o_RHS, T);
  rhsPmlVolume(mesh.NpmlElements, mesh.o_pmlElements, mesh.o_pmlIds,
               o_Q, o_pmlQ, o_RHS, o_pmlRHS, T);

  // compute relaxation terms
  rhsRelaxation(mesh.NnonPmlElements, mesh.o_nonPmlElements, o_Q, o_RHS);
  rhsPmlRelaxation(mesh.NpmlElements, mesh.o_pmlElements, mesh.o_pmlIds,
                   o_Q, o_pmlQ, o_RHS, o_pmlRHS);

  // complete trace halo exchange
  traceHalo->ExchangeFinish(o_Q, 1, ogs_dfloat);

  // compute surface contribution to bns RHS
  rhsSurface(mesh.NnonPmlElements, mesh.o_nonPmlElements, o_Q, o_RHS, T);
  rhsPmlSurface(mesh.NpmlElements, mesh.o_pmlElements, mesh.o_pmlIds,
                o_Q, o_pmlQ, o_RHS, o_pmlRHS, T);
}


//evaluate ODE rhs = f(q,t)
void bns_t::rhsf_MR_pml(occa::memory& o_Q, occa::memory& o_pmlQ,
                        occa::memory& o_RHS, occa::memory& o_pmlRHS,
                        occa::memory& o_fQM, const dfloat T, const int lev){

  // extract q trace halo and start exchange
  multirateTraceHalo[lev]->ExchangeStart(o_fQM, 1, ogs_dfloat);

  // compute volume contribution to bns RHS
  rhsVolume(mesh.mrNnonPmlElements[lev], mesh.o_mrNonPmlElements[lev], o_Q, o_RHS, T);
  rhsPmlVolume(mesh.mrNpmlElements[lev], mesh.o_mrPmlElements[lev], mesh.o_mrPmlIds[lev],
               o_Q, o_pmlQ, o_RHS, o_pmlRHS, T);

  // compute relaxation terms
  rhsRelaxation(mesh.mrNnonPmlElements[lev], mesh.o_mrNonPmlElements[lev], o_Q, o_RHS);
  rhsPmlRelaxation(mesh.mrNpmlElements[lev], mesh.o_mrPmlElements[lev], mesh.o_mrPmlIds[lev],
                   o_Q, o_pmlQ, o_RHS, o_pmlRHS);

  // complete trace halo exchange
  multirateTraceHalo[lev]->ExchangeFinish(o_fQM, 1, ogs_dfloat);

  // compute surface contribution to bns RHS
  rhsSurfaceMR(mesh.mrNnonPmlElements[lev], mesh.o_mrNonPmlElements[lev], o_Q, o_RHS, o_fQM, T);
  rhsPmlSurfaceMR(mesh.mrNpmlElements[lev], mesh.o_mrPmlElements[lev], mesh.o_mrPmlIds[lev],
                  o_Q, o_pmlQ, o_RHS, o_pmlRHS, o_fQM, T);
}

void bns_t::rhsVolume(dlong N, occa::memory& o_ids,
                      occa::memory& o_Q, occa::memory& o_RHS, const dfloat T){

  // compute volume contribution to bns RHS
  if (N)
    volumeKernel(N,
                 o_ids,
                 mesh.o_vgeo,
                 mesh.o_Dmatrices,
                 mesh.o_x,
                 mesh.o_y,
                 mesh.o_z,
                 T,
                 c,
                 nu,
                 o_Q,
                 o_RHS);
}

void bns_t::rhsPmlVolume(dlong N, occa::memory& o_ids, occa::memory& o_pmlids,
                         occa::memory& o_Q, occa::memory& o_pmlQ,
                         occa::memory& o_RHS, occa::memory& o_pmlRHS, const dfloat T){

  // compute volume contribution to bns RHS
  if (N) {
    if (pmlcubature)
      pmlVolumeKernel(N,
                     o_ids,
                     o_pmlids,
                     mesh.o_vgeo,
                     mesh.o_Dmatrices,
                     mesh.o_x,
                     mesh.o_y,
                     mesh.o_z,
                     T,
                     c,
                     nu,
                     o_Q,
                     o_RHS,
                     o_pmlRHS);
    else
      pmlVolumeKernel(N,
                     o_ids,
                     o_pmlids,
                     mesh.o_vgeo,
                     mesh.o_Dmatrices,
                     mesh.o_x,
                     mesh.o_y,
                     mesh.o_z,
                     o_pmlSigma,
                     pmlAlpha,
                     T,
                     c,
                     nu,
                     o_Q,
                     o_pmlQ,
                     o_RHS,
                     o_pmlRHS);
  }
}

void bns_t::rhsRelaxation(dlong N, occa::memory& o_ids,
                          occa::memory& o_Q, occa::memory& o_RHS){

  // compute volume contribution to bns RHS
  if (N)
    relaxationKernel(N,
                     o_ids,
                     mesh.o_vgeo,
                     mesh.o_cubvgeo,
                     mesh.o_cubInterpT,
                     mesh.o_cubProjectT,
                     semiAnalytic,
                     tauInv,
                     o_Q,
                     o_RHS);
}

void bns_t::rhsPmlRelaxation(dlong N, occa::memory& o_ids, occa::memory& o_pmlids,
                             occa::memory& o_Q, occa::memory& o_pmlQ,
                             occa::memory& o_RHS, occa::memory& o_pmlRHS){

  // compute volume contribution to bns RHS
  if (N) {
    if (pmlcubature)
      pmlRelaxationKernel(N,
                         o_ids,
                         o_pmlids,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubInterpT,
                         mesh.o_cubProjectT,
                         o_pmlSigma,
                         pmlAlpha,
                         semiAnalytic,
                         tauInv,
                         o_Q,
                         o_pmlQ,
                         o_RHS,
                         o_pmlRHS);
    else
      pmlRelaxationKernel(N,
                         o_ids,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubInterpT,
                         mesh.o_cubProjectT,
                         semiAnalytic,
                         tauInv,
                         o_Q,
                         o_RHS);
  }
}

void bns_t::rhsSurface(dlong N, occa::memory& o_ids,
                      occa::memory& o_Q, occa::memory& o_RHS, const dfloat T){

  // compute volume contribution to bns RHS
  if (N)
    surfaceKernel(N,
                  o_ids,
                  mesh.o_sgeo,
                  mesh.o_LIFTT,
                  mesh.o_vmapM,
                  mesh.o_vmapP,
                  mesh.o_EToB,
                  mesh.o_x,
                  mesh.o_y,
                  mesh.o_z,
                  T,
                  c,
                  nu,
                  o_Q,
                  o_RHS);
}

void bns_t::rhsPmlSurface(dlong N, occa::memory& o_ids, occa::memory& o_pmlids,
                         occa::memory& o_Q, occa::memory& o_pmlQ,
                         occa::memory& o_RHS, occa::memory& o_pmlRHS, const dfloat T){

  // compute volume contribution to bns RHS
  if (N)
    pmlSurfaceKernel(N,
                    o_ids,
                    o_pmlids,
                    mesh.o_sgeo,
                    mesh.o_LIFTT,
                    mesh.o_vmapM,
                    mesh.o_vmapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    T,
                    c,
                    nu,
                    o_Q,
                    o_RHS,
                    o_pmlRHS);
}

void bns_t::rhsSurfaceMR(dlong N, occa::memory& o_ids,
                         occa::memory& o_Q, occa::memory& o_RHS,
                         occa::memory& o_fQM, const dfloat T){

  // compute volume contribution to bns RHS
  if (N)
    surfaceKernel(N,
                  o_ids,
                  mesh.o_sgeo,
                  mesh.o_LIFTT,
                  mesh.o_vmapM,
                  mesh.o_mapP,
                  mesh.o_EToB,
                  mesh.o_x,
                  mesh.o_y,
                  mesh.o_z,
                  T,
                  c,
                  nu,
                  o_Q,
                  o_fQM,
                  o_RHS);
}

void bns_t::rhsPmlSurfaceMR(dlong N, occa::memory& o_ids, occa::memory& o_pmlids,
                         occa::memory& o_Q, occa::memory& o_pmlQ,
                         occa::memory& o_RHS, occa::memory& o_pmlRHS,
                         occa::memory& o_fQM, const dfloat T){

  // compute volume contribution to bns RHS
  if (N)
    pmlSurfaceKernel(N,
                    o_ids,
                    o_pmlids,
                    mesh.o_sgeo,
                    mesh.o_LIFTT,
                    mesh.o_vmapM,
                    mesh.o_mapP,
                    mesh.o_EToB,
                    mesh.o_x,
                    mesh.o_y,
                    mesh.o_z,
                    T,
                    c,
                    nu,
                    o_Q,
                    o_fQM,
                    o_RHS,
                    o_pmlRHS);
}
