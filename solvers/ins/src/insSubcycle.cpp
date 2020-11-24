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

#include "ins.hpp"

subcycler_t::subcycler_t(ins_t& ins):
  solver_t(ins.platform, ins.settings), mesh(ins.mesh) {

  NVfields = ins.NVfields;
  nu = ins.nu;
  cubature = ins.cubature;
  vTraceHalo = ins.vTraceHalo;
  advectionVolumeKernel = ins.advectionVolumeKernel;
  advectionSurfaceKernel = ins.advectionSurfaceKernel;
  advectionInterpolationKernel = ins.advectionInterpolationKernel;

  int mOrder = 6; // HACK
  o_cUh = platform.malloc(mOrder*ins.mesh.cubNp*ins.mesh.Nelements*NVfields*sizeof(dfloat));
  //  o_GUh = platform.malloc(mOrder*ins.mesh.Np*ins.mesh.Nelements*NVfields*sizeof(dfloat));
}

//evaluate ODE rhs = f(q,t)
void subcycler_t::rhsf(occa::memory& o_U, occa::memory& o_RHS, const dfloat T){

  dfloat c0 = 1, c1 = 0, c2 = 0;
  const dfloat t0 = T;
  const dfloat t1 = T-dt;
  const dfloat t2 = T-2*dt;
  
  switch(order){
  case 0:
    c0 = 1; c1 = 0; c2 = 0; break;
  case 1:
    c0 = (T-t1)/(t0-t1);  c1 = (T-t0)/(t1-t0);  c2 = 0;  break;
  case 2:
    c0 = (T-t1)*(T-t2)/((t0-t1)*(t0-t2));
    c1 = (T-t0)*(T-t2)/((t1-t0)*(t1-t2));
    c2 = (T-t0)*(T-t1)/((t2-t0)*(t2-t1));
    break;
  default:
    printf("WEIRD ORDER: %d\n", order);
    exit(-1);
  }

  dlong fieldOffset = mesh.Nelements*mesh.Np*NVfields;
  dlong offset0 = ((shiftIndex+0)%maxOrder)*fieldOffset;
  dlong offset1 = ((shiftIndex+1)%maxOrder)*fieldOffset;
  dlong offset2 = ((shiftIndex+2)%maxOrder)*fieldOffset;

  // TW: need to time interpolate GUh to GUe instead of Uh to Ue
  
  //interpolate velocity history for advective field (halo elements first)
  if(mesh.NhaloElements){
    subCycleAdvectionKernel(mesh.NhaloElements,
			    mesh.o_haloElementIds,
			    mesh.Np,
			    order,
			    offset0, offset1, offset2,
			    c0, c1, c2,
			    o_Uh,
			    o_Ue);
  }
  
  // extract Ue halo
  vTraceHalo->ExchangeStart(o_Ue, 1, ogs_dfloat);
  
  if(mesh.NinternalElements){
    subCycleAdvectionKernel(mesh.NinternalElements,
			    mesh.o_internalElementIds,
			    mesh.Np,
			    order,
			    offset0, offset1, offset2,
			    c0, c1, c2,
			    o_Uh,
			    o_Ue);

    if(cubature){
      dlong cubFieldOffset = mesh.Nelements*mesh.cubNp*NVfields;
      dlong cubOffset0 = ((shiftIndex+0)%maxOrder)*cubFieldOffset;
      dlong cubOffset1 = ((shiftIndex+1)%maxOrder)*cubFieldOffset;
      dlong cubOffset2 = ((shiftIndex+2)%maxOrder)*cubFieldOffset;
      
      subCycleAdvectionKernel(mesh.NinternalElements,
			      mesh.o_internalElementIds,
			      mesh.cubNp,
			      order,
			      cubOffset0, cubOffset1, cubOffset2,
			      c0, c1, c2,
			      o_cUh,
			      o_cUe);
    }
  }


  // finish exchange of Ue
  vTraceHalo->ExchangeFinish(o_Ue, 1, ogs_dfloat);

  // extract u halo on DEVICE
  vTraceHalo->ExchangeStart(o_U, 1, ogs_dfloat);

  if (cubature)
    advectionVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_cubvgeo,
                         mesh.o_cubD,
                         mesh.o_cubPDT,
                         mesh.o_cubInterp,
                         mesh.o_cubProject,
                         o_cUe,
                         o_U,
                         o_RHS);
  else
    advectionVolumeKernel(mesh.Nelements,
                         mesh.o_vgeo,
                         mesh.o_D,
                         o_Ue,
                         o_U,
                         o_RHS);

  vTraceHalo->ExchangeFinish(o_U, 1, ogs_dfloat);

  if (cubature)
    advectionSurfaceKernel(mesh.Nelements,
                          mesh.o_vgeo,
                          mesh.o_cubsgeo,
                          mesh.o_intInterp,
                          mesh.o_intLIFT,
                          mesh.o_vmapM,
                          mesh.o_vmapP,
                          mesh.o_EToB,
                          T,
                          mesh.o_intx,
                          mesh.o_inty,
                          mesh.o_intz,
                          nu,
                          o_Ue,
                          o_U,
                          o_RHS);
  else
    advectionSurfaceKernel(mesh.Nelements,
                          mesh.o_sgeo,
                          mesh.o_LIFT,
                          mesh.o_vmapM,
                          mesh.o_vmapP,
                          mesh.o_EToB,
                          T,
                          mesh.o_x,
                          mesh.o_y,
                          mesh.o_z,
                          nu,
                          o_Ue,
                          o_U,
                          o_RHS);
}
